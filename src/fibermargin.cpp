#include <Rcpp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <numeric>
#include <thread>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace fiber_margin {

constexpr double kGoldenRatioConjugate = 0.6180339887498949;
constexpr double kSqrtTwoConjugate = 0.4142135623730951;
constexpr double kSqrtThreeConjugate = 0.7320508075688772;
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr int kViews2D = 9;
constexpr int kViews3D = 12;
constexpr double kStandardTransportRange = 5.0;
constexpr double kShortTransportRange = 2.5;
constexpr int kMaxWorkers = 64;

struct Settings {
  int threads = 1;
  double transport_range = kStandardTransportRange;
};

int list_int(const List& control, const char* name, int fallback) {
  if (!control.containsElementNamed(name)) return fallback;
  return as<int>(control[name]);
}

Settings parse_settings(const List& control) {
  const int known = control.containsElementNamed("threads") +
                    control.containsElementNamed("short_transport");
  if (control.size() != known) {
    stop("FiberMargin has no tuning controls beyond the documented presets.");
  }
  Settings settings;
  settings.threads = list_int(control, "threads", settings.threads);
  if (settings.threads < 1 || settings.threads > kMaxWorkers) {
    stop("FiberMargin threads must be between 1 and 64.");
  }
  const int short_transport = list_int(control, "short_transport", 0);
  if (short_transport != 0 && short_transport != 1) stop("Invalid transport preset.");
  if (short_transport) settings.transport_range = kShortTransportRange;
  return settings;
}

double quantile_linear(std::vector<double> values, double probability) {
  if (values.empty()) return 0.0;
  const double position = probability * static_cast<double>(values.size() - 1);
  const std::size_t lower = static_cast<std::size_t>(std::floor(position));
  const std::size_t upper = static_cast<std::size_t>(std::ceil(position));
  const double fraction = position - static_cast<double>(lower);
  const auto pivot = values.begin() + upper;
  std::nth_element(values.begin(), pivot, values.end());
  const double lower_value = lower == upper ? *pivot : *std::max_element(values.begin(), pivot);
  return lower_value * (1.0 - fraction) + *pivot * fraction;
}

double median_inplace(std::vector<double>& values) {
  if (values.empty()) return 0.0;
  const std::size_t middle = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + middle, values.end());
  const double upper = values[middle];
  if (values.size() % 2 == 1) return upper;
  const double lower = *std::max_element(values.begin(), values.begin() + middle);
  return 0.5 * (lower + upper);
}

struct CoordinateBlock {
  std::vector<double> values;
  int dimensions;
};

CoordinateBlock robust_unit_coordinates(
    const NumericMatrix& xy, const std::vector<int>& rows) {
  std::vector<int> active_axes;
  for (int axis = 0; axis < xy.ncol(); ++axis) {
    double low = xy(rows.front(), axis);
    double high = low;
    for (std::size_t i = 1; i < rows.size(); ++i) {
      low = std::min(low, xy(rows[i], axis));
      high = std::max(high, xy(rows[i], axis));
    }
    if (low < high) active_axes.push_back(axis);
  }
  const int dimensions = static_cast<int>(active_axes.size());
  if (dimensions < 2) {
    stop("Each sample needs at least two varying coordinate axes.");
  }
  const int size = static_cast<int>(rows.size());
  CoordinateBlock output;
  output.dimensions = dimensions;
  output.values.resize(static_cast<std::size_t>(size) * dimensions);
  for (int dimension = 0; dimension < dimensions; ++dimension) {
    const int axis = active_axes[dimension];
    std::vector<double> values(size);
    for (int i = 0; i < size; ++i) values[i] = xy(rows[i], axis);
    const auto extent = std::minmax_element(values.begin(), values.end());
    double axis_low = *extent.first;
    double axis_high = *extent.second;
    const double magnitude = std::max(std::abs(axis_low), std::abs(axis_high));
    // Shift extreme exponents before forming differences or interpolating.
    if (magnitude > std::numeric_limits<double>::max() / 4 ||
        magnitude < 4 * std::numeric_limits<double>::min()) {
      int exponent = 0;
      std::frexp(magnitude, &exponent);
      for (double& value : values) value = std::scalbn(value, -exponent);
      axis_low = std::scalbn(axis_low, -exponent);
      axis_high = std::scalbn(axis_high, -exponent);
    }
    double low = quantile_linear(values, 0.01);
    double high = quantile_linear(values, 0.99);
    if (high <= low) {
      low = axis_low;
      high = axis_high;
    }
    const double span = high - low;
    for (int i = 0; i < size; ++i) {
      const double value = (values[i] - low) / span;
      output.values[static_cast<std::size_t>(i) * dimensions + dimension] =
        std::max(-0.1, std::min(1.1, value));
    }
  }
  return output;
}

std::vector<double> rotation_matrix(int dimensions, int view) {
  if (dimensions == 2) {
    const double angle = kPi * std::fmod(view * kGoldenRatioConjugate, 1.0);
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    return {cosine, -sine, sine, cosine};
  }

  const double u1 = std::fmod(0.5 + view * kSqrtTwoConjugate, 1.0);
  const double u2 = std::fmod(0.25 + view * kGoldenRatioConjugate, 1.0);
  const double u3 = std::fmod(0.75 + view * kSqrtThreeConjugate, 1.0);
  const double x = std::sqrt(1.0 - u1) * std::sin(2.0 * kPi * u2);
  const double y = std::sqrt(1.0 - u1) * std::cos(2.0 * kPi * u2);
  const double z = std::sqrt(u1) * std::sin(2.0 * kPi * u3);
  const double w = std::sqrt(u1) * std::cos(2.0 * kPi * u3);
  return {
    1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w),
    2.0 * (x * z + y * w), 2.0 * (x * y + z * w),
    1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - x * w),
    2.0 * (x * z - y * w), 2.0 * (y * z + x * w),
    1.0 - 2.0 * (x * x + y * y)
  };
}

std::vector<double> path_unit(const std::vector<double>& coordinates,
                              int size, int dimensions, int view) {
  const std::vector<double> rotation = rotation_matrix(dimensions, view);
  std::vector<double> transformed(static_cast<std::size_t>(size) * dimensions);
  std::vector<double> low(dimensions, std::numeric_limits<double>::infinity());
  std::vector<double> high(dimensions, -std::numeric_limits<double>::infinity());
  for (int i = 0; i < size; ++i) {
    for (int destination = 0; destination < dimensions; ++destination) {
      double value = 0.0;
      for (int source = 0; source < dimensions; ++source) {
        value += coordinates[static_cast<std::size_t>(i) * dimensions + source] *
          rotation[static_cast<std::size_t>(source) * dimensions + destination];
      }
      transformed[static_cast<std::size_t>(i) * dimensions + destination] = value;
      low[destination] = std::min(low[destination], value);
      high[destination] = std::max(high[destination], value);
    }
  }
  for (int i = 0; i < size; ++i) {
    for (int dimension = 0; dimension < dimensions; ++dimension) {
      const std::size_t index = static_cast<std::size_t>(i) * dimensions + dimension;
      transformed[index] = std::max(0.0, std::min(
        1.0, (transformed[index] - low[dimension]) /
          std::max(high[dimension] - low[dimension], 1e-12)));
    }
  }
  return transformed;
}

// Spread the 16 low bits into alternating bit lanes for the Hilbert transpose.
template <int dimensions>
std::uint64_t spread_bits(std::uint64_t value) {
  if constexpr (dimensions == 2) {
    value = (value | value << 8) & 0x00ff00ffULL;
    value = (value | value << 4) & 0x0f0f0f0fULL;
    value = (value | value << 2) & 0x33333333ULL;
    return (value | value << 1) & 0x55555555ULL;
  } else {
    value = (value | value << 16) & 0x001f0000ff0000ffULL;
    value = (value | value << 8) & 0x100f00f00f00f00fULL;
    value = (value | value << 4) & 0x10c30c30c30c30c3ULL;
    return (value | value << 2) & 0x1249249249249249ULL;
  }
}

template <int dimensions>
std::uint64_t hilbert_code(const double* coordinate) {
  std::array<std::uint64_t, dimensions> transformed;
  constexpr int bits = 16;
  constexpr std::uint64_t maximum = (1ULL << bits) - 1;
  for (int dimension = 0; dimension < dimensions; ++dimension) {
    const double clipped = std::max(0.0, std::min(1.0 - 1e-12, coordinate[dimension]));
    transformed[dimension] =
      static_cast<std::uint64_t>(std::floor(clipped * maximum));
  }

  for (std::uint64_t q = 1ULL << (bits - 1); q > 1; q >>= 1) {
    const std::uint64_t p = q - 1;
    for (int axis = 0; axis < dimensions; ++axis) {
      const std::uint64_t invert = 0ULL - ((transformed[axis] & q) != 0);
      const std::uint64_t exchange = (transformed[0] ^ transformed[axis]) & p & ~invert;
      transformed[0] ^= (invert & p) | exchange;
      transformed[axis] ^= exchange;
    }
  }
  for (int axis = 1; axis < dimensions; ++axis) {
    transformed[axis] ^= transformed[axis - 1];
  }
  std::uint64_t correction = transformed[dimensions - 1] >> 1;
  correction ^= correction >> 1;
  correction ^= correction >> 2;
  correction ^= correction >> 4;
  correction ^= correction >> 8;
  std::uint64_t code = 0;
  for (int axis = 0; axis < dimensions; ++axis) {
    code |= spread_bits<dimensions>(transformed[axis] ^ correction) << (dimensions - 1 - axis);
  }
  return code;
}

struct Geometry {
  std::vector<int> order;
  std::vector<std::pair<int, int>> coincident;
  std::vector<double> isolation;
  std::vector<double> decay;
};

int compare_coordinates(const std::vector<double>& coordinates, int dimensions,
                        int left, int right) {
  const auto* a = coordinates.data() + static_cast<std::size_t>(left) * dimensions;
  const auto* b = coordinates.data() + static_cast<std::size_t>(right) * dimensions;
  for (int axis = 0; axis < dimensions; ++axis) {
    if (a[axis] < b[axis]) return -1;
    if (a[axis] > b[axis]) return 1;
  }
  return 0;
}

template <int dimensions>
Geometry path_geometry(const std::vector<double>& coordinates, int size, int view,
                       double transport_range) {
  Geometry geometry;
  geometry.order.resize(size);
  geometry.isolation.assign(size, 1.0);
  std::iota(geometry.order.begin(), geometry.order.end(), 0);

  std::vector<std::uint64_t> code(size);
  {
    const std::vector<double> unit = path_unit(coordinates, size, dimensions, view);
    for (int i = 0; i < size; ++i) {
      code[i] = hilbert_code<dimensions>(&unit[static_cast<std::size_t>(i) * dimensions]);
    }
  }
  std::sort(geometry.order.begin(), geometry.order.end(), [&](int left, int right) {
    if (code[left] != code[right]) return code[left] < code[right];
    return compare_coordinates(coordinates, dimensions, left, right) < 0;
  });

  if (size < 2) return geometry;
  std::vector<double> step(size - 1);
  for (int position = 1; position < size; ++position) {
    const int left = geometry.order[position - 1];
    const int right = geometry.order[position];
    double squared = 0.0;
    for (int dimension = 0; dimension < dimensions; ++dimension) {
      const double delta =
        coordinates[static_cast<std::size_t>(left) * dimensions + dimension] -
        coordinates[static_cast<std::size_t>(right) * dimensions + dimension];
      squared += delta * delta;
    }
    step[position - 1] = std::sqrt(squared);
  }
  std::vector<double> positive;
  positive.reserve(step.size());
  for (double value : step) if (value > 1e-12) positive.push_back(value);
  const double reference = positive.empty() ? 1.0 : median_inplace(positive);
  const auto local_scale = [&](int edge) { return std::max(step[edge], 1e-12); };
  for (int position = 0; position < size; ++position) {
    const double scale = position == 0 ? local_scale(0) : position == size - 1 ?
      local_scale(size - 2) : std::sqrt(local_scale(position - 1) * local_scale(position));
    const double gap = scale / reference;
    geometry.isolation[geometry.order[position]] = std::max(1.0, gap);
  }
  for (int begin = 0; begin < size - 1;) {
    int end = begin + 1;
    while (end < size && step[end - 1] == 0.0 &&
           compare_coordinates(coordinates, dimensions,
                               geometry.order[begin], geometry.order[end]) == 0) ++end;
    if (end > begin + 1) {
      geometry.coincident.emplace_back(begin, end);
      const double left = begin > 0 ? local_scale(begin - 1) :
        (end < size ? local_scale(end - 1) : reference);
      const double right = end < size ? local_scale(end - 1) : left;
      const double gap = std::max(1.0, std::sqrt(left * right) / reference);
      for (int i = begin; i < end; ++i) geometry.isolation[geometry.order[i]] = gap;
    }
    begin = end;
  }
  for (double& edge : step) edge = std::exp(-(edge / reference) / transport_range);
  geometry.decay = std::move(step);
  return geometry;
}

// A class is spatially enclosed at a path position when its transported mass
// arrives from both directions. The squared sum of directional amplitudes
// keeps unilateral support while giving an explicit bonus to balanced support:
// (sqrt(left) + sqrt(right))^2.
void two_sided_enclosure_field(const std::vector<int>& labels, int classes,
                               const Geometry& geometry, std::vector<double>& field) {
  const int size = static_cast<int>(labels.size());
  const std::size_t width = static_cast<std::size_t>(size) * classes;
  field.resize(width);
  std::vector<double> state(classes, 0.0);
  const auto decay_at = [&](int index, bool reverse) {
    if (reverse) return index == size - 1 ? 1.0 : geometry.decay[index];
    return index == 0 ? 1.0 : geometry.decay[index - 1];
  };
  const auto sweep = [&](bool reverse_direction, const auto& visit) {
    if (geometry.coincident.empty()) {
      for (int pos = 0; pos < size; ++pos) {
        const int index = reverse_direction ? size - 1 - pos : pos;
        const int row = geometry.order[index];
        const double decay = decay_at(index, reverse_direction);
        const std::size_t offset = static_cast<std::size_t>(row) * classes;
        for (int cls = 0; cls < classes; ++cls) {
          state[cls] *= decay;
          visit(offset + cls, state[cls]);
        }
        state[labels[row]] += 1.0;
      }
      return;
    }
    std::vector<int> counts(classes);
    int group = reverse_direction ? static_cast<int>(geometry.coincident.size()) - 1 : 0;
    for (int pos = 0; pos < size;) {
      const int index = reverse_direction ? size - 1 - pos : pos;
      const int row = geometry.order[index];
      const double decay = decay_at(index, reverse_direction);
      int length = 1;
      if (group >= 0 && group < static_cast<int>(geometry.coincident.size())) {
        const auto& block = geometry.coincident[group];
        if (index == (reverse_direction ? block.second - 1 : block.first)) {
          length = block.second - block.first;
          group += reverse_direction ? -1 : 1;
        }
      }
      const auto row_at = [&](int j) {
        return geometry.order[index + (reverse_direction ? -j : j)];
      };
      if (length == 1) {
        const std::size_t offset = static_cast<std::size_t>(row) * classes;
        for (int cls = 0; cls < classes; ++cls) {
          state[cls] *= decay;
          visit(offset + cls, state[cls]);
        }
        state[labels[row]] += 1.0;
        ++pos;
        continue;
      }
      std::fill(counts.begin(), counts.end(), 0);
      for (int j = 0; j < length; ++j) ++counts[labels[row_at(j)]];
      for (int cls = 0; cls < classes; ++cls) state[cls] *= decay;
      // Split peers' mass equally between directions, excluding the query itself.
      for (int j = 0; j < length; ++j) {
        const int peer = row_at(j);
        const int observed = labels[peer];
        const std::size_t offset = static_cast<std::size_t>(peer) * classes;
        for (int cls = 0; cls < classes; ++cls) {
          visit(offset + cls, state[cls] + 0.5 * (counts[cls] - (cls == observed)));
        }
      }
      for (int cls = 0; cls < classes; ++cls) state[cls] += counts[cls];
      pos += length;
    }
  };
  sweep(false, [&](std::size_t index, double left) { field[index] = left; });
  std::fill(state.begin(), state.end(), 0.0);
  sweep(true, [&](std::size_t index, double reverse) {
    const double forward = field[index];
    field[index] = forward + reverse + 2.0 * std::sqrt(forward * reverse);
  });
}

bool has_spatial_extent(const std::vector<double>& coordinates, int size,
                        int dimensions) {
  if (size < 2) return false;
  for (int dimension = 0; dimension < dimensions; ++dimension) {
    const double first = coordinates[dimension];
    for (int row = 1; row < size; ++row) {
      if (coordinates[static_cast<std::size_t>(row) * dimensions + dimension] != first) {
        return true;
      }
    }
  }
  return false;
}

struct SampleResult {
  std::vector<int> labels;
  std::vector<int> candidate;
  std::vector<double> margin_score;
  std::vector<double> required;
  std::vector<double> dispersion;
  std::vector<double> isolation;
};

// Each spatial route supplies pointwise candidate evidence. The final label
// decision remains separate from its construction, so binary ballots and
// multiclass path transport share one selective margin rule.
struct MarginEvidence {
  std::vector<int> candidate;
  std::vector<double> margin_score;
  std::vector<double> required;
  std::vector<double> dispersion;
  std::vector<double> isolation;
};

struct NeighborNode {
  int row;
  int axis;
  int left;
  int right;
};

struct ThreadLauncher {
  template <typename Task>
  void operator()(std::vector<std::thread>& pool, Task task) const {
    pool.emplace_back(std::move(task));
  }
};

template <typename Task, typename Launch = ThreadLauncher>
void run_workers(int workers, const Task& task, const Launch& launch = Launch{}) {
  if (workers == 1) {
    task(0);
    return;
  }
  std::vector<std::exception_ptr> errors(workers);
  std::vector<std::thread> pool;
  pool.reserve(workers);
  struct JoinOnExit {
    std::vector<std::thread>& pool;
    ~JoinOnExit() {
      for (auto& worker : pool) if (worker.joinable()) worker.join();
    }
  };
  {
    JoinOnExit join{pool};
    for (int i = 0; i < workers; ++i) {
      launch(pool, [&, i]() {
        try {
          task(i);
        } catch (...) {
          errors[i] = std::current_exception();
        }
      });
    }
  }
  for (const auto& error : errors) if (error) std::rethrow_exception(error);
}

class NeighborTree {
 public:
  NeighborTree(const std::vector<double>& coordinates, int dimensions,
               const std::vector<int>& rows, const std::vector<int>& labels)
      : coordinates_(coordinates), dimensions_(dimensions), order_(rows), rank_(rows.size()) {
    std::sort(order_.begin(), order_.end(), [&](int left, int right) {
      const int comparison = compare_coordinates(coordinates, dimensions, left, right);
      return comparison == 0 ? labels[left] < labels[right] : comparison < 0;
    });
    by_rank_ = order_;
    for (std::size_t j = 0; j < order_.size(); ++j) rank_[order_[j]] = static_cast<int>(j);
    nodes_.reserve(rows.size());
    root_ = build(0, static_cast<int>(order_.size()), 0);
  }

  void binary_counts(int target, int count, const std::vector<int>& labels,
                     std::vector<std::pair<double, int> >& heap,
                     int& count_zero, int& count_one) const {
    heap.clear();
    search(root_, target, count, heap);
    count_zero = 0;
    count_one = 0;
    for (const std::pair<double, int>& entry : heap) {
      if (labels[by_rank_[entry.second]] == 0) {
        ++count_zero;
      } else {
        ++count_one;
      }
    }
  }

 private:
  const std::vector<double>& coordinates_;
  int dimensions_;
  int root_;
  std::vector<int> order_;
  std::vector<int> rank_;
  std::vector<int> by_rank_;
  std::vector<NeighborNode> nodes_;

  double coordinate(int row, int axis) const {
    return coordinates_[static_cast<std::size_t>(row) * dimensions_ + axis];
  }

  int build(int begin, int end, int depth) {
    if (begin >= end) return -1;
    const int axis = depth % dimensions_;
    const int middle = begin + (end - begin) / 2;
    std::nth_element(
      order_.begin() + begin, order_.begin() + middle, order_.begin() + end,
      [&](int left, int right) { return coordinate(left, axis) < coordinate(right, axis); }
    );
    const int node = static_cast<int>(nodes_.size());
    nodes_.push_back({order_[middle], axis, -1, -1});
    nodes_[node].left = build(begin, middle, depth + 1);
    nodes_[node].right = build(middle + 1, end, depth + 1);
    return node;
  }

  double squared_distance(int left, int right) const {
    double total = 0.0;
    for (int axis = 0; axis < dimensions_; ++axis) {
      const double delta = coordinate(left, axis) - coordinate(right, axis);
      total += delta * delta;
    }
    return total;
  }

  void search(int node_index, int target, int count,
              std::vector<std::pair<double, int> >& heap) const {
    if (node_index < 0) return;
    const NeighborNode& node = nodes_[node_index];
    const double difference = coordinate(target, node.axis) - coordinate(node.row, node.axis);
    const int near_child = difference <= 0.0 ? node.left : node.right;
    const int far_child = difference <= 0.0 ? node.right : node.left;
    search(near_child, target, count, heap);
    if (node.row != target) {
      const double distance = squared_distance(target, node.row);
      const auto entry = std::make_pair(distance, rank_[node.row]);
      if (static_cast<int>(heap.size()) < count) {
        heap.push_back(entry);
        std::push_heap(heap.begin(), heap.end());
      } else if (entry < heap.front()) {
        std::pop_heap(heap.begin(), heap.end());
        heap.back() = entry;
        std::push_heap(heap.begin(), heap.end());
      }
    }
    const double bound = static_cast<int>(heap.size()) < count ?
      std::numeric_limits<double>::infinity() : heap.front().first;
    if (difference * difference <= bound) search(far_child, target, count, heap);
  }
};

template <typename Function>
void parallel_ranges(int size, int threads, const Function& function) {
  const int workers = std::max(1, std::min(size, threads));
  run_workers(workers, [&](int worker) {
    const int begin = static_cast<std::int64_t>(size) * worker / workers;
    const int end = static_cast<std::int64_t>(size) * (worker + 1) / workers;
    function(begin, end);
  });
}

SampleResult apply_margin(const std::vector<int>& labels,
                          MarginEvidence evidence,
                          const Settings& settings) {
  const int size = static_cast<int>(labels.size());
  SampleResult result;
  result.labels = labels;
  result.candidate = std::move(evidence.candidate);
  result.margin_score = std::move(evidence.margin_score);
  result.required = std::move(evidence.required);
  result.dispersion = std::move(evidence.dispersion);
  result.isolation = std::move(evidence.isolation);
  // A short admission pass costs less than starting a new thread pool.
  const int workers = size < 4096 ? 1 : settings.threads;
  parallel_ranges(size, workers, [&](int begin, int end) {
    for (int row = begin; row < end; ++row) {
      if (result.candidate[row] != labels[row] &&
          result.margin_score[row] >= result.required[row]) {
        result.labels[row] = result.candidate[row];
      }
    }
  });
  return result;
}

MarginEvidence binary_ballot_evidence(
    const std::vector<double>& coordinates, const std::vector<int>& labels,
    int dimensions, const Settings& settings) {
  const int size = static_cast<int>(labels.size());
  MarginEvidence evidence;
  evidence.candidate = labels;
  evidence.margin_score.assign(size, 0.0);
  // Integer ballot scores require one net local-volume vote to displace the
  // supplied category. Exact ties therefore retain that category.
  evidence.required.assign(size, 1.0);
  evidence.dispersion.assign(size, 0.0);
  evidence.isolation.assign(size, 1.0);
  if (size < 3) return evidence;

  std::vector<int> rows(size);
  std::iota(rows.begin(), rows.end(), 0);
  const NeighborTree tree(coordinates, dimensions, rows, labels);
  const int local_mass = std::min(size - 1, 1 << (dimensions + 4));
  parallel_ranges(size, settings.threads, [&](int begin, int end) {
    std::vector<std::pair<double, int> > heap;
    heap.reserve(local_mass);
    for (int row = begin; row < end; ++row) {
      int count_zero = 0;
      int count_one = 0;
      tree.binary_counts(row, local_mass, labels, heap, count_zero, count_one);
      // A tied local volume has no preferred rival. Retaining the observed
      // class makes the binary evidence obey the same candidate-first rule as
      // the multiclass branch and prevents diagnostic-only tie choices.
      const int rival = count_zero > count_one ? 0 :
        (count_one > count_zero ? 1 : labels[row]);
      const int rival_count = rival == 0 ? count_zero : count_one;
      const int observed_count = labels[row] == 0 ? count_zero : count_one;
      evidence.candidate[row] = rival;
      evidence.margin_score[row] = rival == labels[row] ? 0.0 :
        static_cast<double>(rival_count - observed_count);
    }
  });
  return evidence;
}

SampleResult refine_sample(const std::vector<double>& coordinates,
                           const std::vector<int>& labels, int dimensions,
                           int classes, const Settings& settings) {
  const int size = static_cast<int>(labels.size());
  MarginEvidence evidence;
  evidence.candidate = labels;
  evidence.margin_score.assign(size, 0.0);
  evidence.required.assign(size, 0.0);
  evidence.dispersion.assign(size, 0.0);
  evidence.isolation.assign(size, 1.0);
  // With no coordinate extent, a path order is arbitrary. Retaining the
  // supplied field avoids creating order-dependent spatial evidence.
  if (!has_spatial_extent(coordinates, size, dimensions)) {
    return apply_margin(labels, std::move(evidence), settings);
  }
  if (size < 20) return apply_margin(labels, std::move(evidence), settings);
  std::vector<double> count(classes, 0.0);
  int observed_classes = 0;
  for (int label : labels) {
    if (count[label] == 0.0) ++observed_classes;
    count[label] += 1.0;
  }
  if (observed_classes < 2) {
    return apply_margin(labels, std::move(evidence), settings);
  }
  if (classes == 2) {
    return apply_margin(
      labels, binary_ballot_evidence(coordinates, labels, dimensions, settings), settings);
  }

  const int charts = dimensions == 2 ? kViews2D : kViews3D;
  const std::size_t row_class = static_cast<std::size_t>(size) * classes;
  std::vector<float> chart_action(static_cast<std::size_t>(charts) * row_class);
  std::vector<float> chart_isolation(static_cast<std::size_t>(charts) * size);

  checkUserInterrupt();
  const int thread_count = std::min(settings.threads, charts);
  run_workers(thread_count, [&](int thread) {
    std::vector<double> field;
    for (int view = thread; view < charts; view += thread_count) {
      const Geometry geometry = dimensions == 2 ?
        path_geometry<2>(coordinates, size, view, settings.transport_range) :
        path_geometry<3>(coordinates, size, view, settings.transport_range);
      two_sided_enclosure_field(labels, classes, geometry, field);
      float* action = chart_action.data() + static_cast<std::size_t>(view) * row_class;
      for (int row = 0; row < size; ++row) {
        chart_isolation[static_cast<std::size_t>(view) * size + row] =
          static_cast<float>(geometry.isolation[row]);
        const std::size_t offset = static_cast<std::size_t>(row) * classes;
        double total = 0.0;
        for (int cls = 0; cls < classes; ++cls) total += field[offset + cls];
        total = std::max(total, 1e-12);
        const double observed = field[offset + labels[row]] / total;
        for (int cls = 0; cls < classes; ++cls) {
          const std::size_t index = offset + cls;
          action[index] = static_cast<float>(field[index] / total - observed);
        }
      }
    }
  });

  std::vector<double> positive;
  for (double value : count) if (value > 0.0) positive.push_back(value);
  const double reference = std::max(median_inplace(positive), 1.0);
  std::vector<double> deficit(classes, 1.0);
  for (int cls = 0; cls < classes; ++cls) {
    deficit[cls] = std::sqrt(reference / std::max(count[cls], 1.0));
  }

  checkUserInterrupt();
  const double sqrt_charts = std::sqrt(static_cast<double>(charts));
  parallel_ranges(size, settings.threads, [&](int begin, int end) {
    std::vector<double> values(charts);
    std::vector<double> mean(classes);
    for (int row = begin; row < end; ++row) {
      for (int atlas = 0; atlas < charts; ++atlas) {
        values[atlas] = chart_isolation[static_cast<std::size_t>(atlas) * size + row];
      }
      evidence.isolation[row] = median_inplace(values);

      const std::size_t offset = static_cast<std::size_t>(row) * classes;
      std::fill(mean.begin(), mean.end(), 0.0);
      for (int atlas = 0; atlas < charts; ++atlas) {
        const float* contrasts = chart_action.data() +
          static_cast<std::size_t>(atlas) * row_class + offset;
        for (int cls = 0; cls < classes; ++cls) mean[cls] += contrasts[cls];
      }
      for (double& value : mean) value /= static_cast<double>(charts);
      int rival = labels[row];
      for (int cls = 0; cls < classes; ++cls) {
        if (mean[cls] > mean[rival]) rival = cls;
      }

      // The candidate uses mean contrasts only; dispersion is needed just for its barrier.
      for (int atlas = 0; atlas < charts; ++atlas) {
        values[atlas] = chart_action[static_cast<std::size_t>(atlas) * row_class + offset + rival];
      }
      double variance = 0.0;
      const double selected_mean = mean[rival];
      for (double value : values) {
        const double delta = value - selected_mean;
        variance += delta * delta;
      }
      variance /= static_cast<double>(charts);
      const double dispersion = std::sqrt(variance) / sqrt_charts;
      evidence.candidate[row] = rival;
      evidence.margin_score[row] = selected_mean;
      evidence.required[row] = dispersion / deficit[rival] * evidence.isolation[row];
      evidence.dispersion[row] = dispersion;
    }
  });
  return apply_margin(labels, std::move(evidence), settings);
}

}  // namespace fiber_margin

extern "C" SEXP _fibermargin_fiber_margin_cpp(SEXP xy_s, SEXP labels_s,
                                              SEXP samples_s, SEXP control_s) {
  BEGIN_RCPP
  NumericMatrix xy(xy_s);
  IntegerVector labels_r(labels_s);
  IntegerVector samples_r(samples_s);
  List control(control_s);
  const int n = xy.nrow();
  const int dimensions = xy.ncol();
  if (dimensions < 2 || dimensions > 3) stop("FiberMargin requires two or three coordinates.");
  if (labels_r.size() != n || samples_r.size() != n) stop("FiberMargin input lengths differ.");

  int classes = 0;
  for (int value : labels_r) {
    if (value == NA_INTEGER || value < 1) stop("FiberMargin labels must be observed factors.");
    classes = std::max(classes, value);
  }
  const fiber_margin::Settings settings = fiber_margin::parse_settings(control);
  IntegerVector output = clone(labels_r);
  IntegerVector candidate = clone(labels_r);
  NumericVector margin_score(n, 0.0);
  NumericVector required(n, 0.0);
  NumericVector dispersion(n, 0.0);
  NumericVector isolation(n, 1.0);

  int max_sample = 0;
  for (int code : samples_r) {
    if (code < 1 || code > n) stop("FiberMargin requires valid specimen codes.");
    max_sample = std::max(max_sample, code);
  }
  std::vector<std::vector<int>> groups(max_sample);
  for (int i = 0; i < n; ++i) groups[samples_r[i] - 1].push_back(i);
  std::vector<int> sample_levels;
  for (int i = 0; i < max_sample; ++i) if (!groups[i].empty()) sample_levels.push_back(i + 1);
  IntegerVector dimensions_used(sample_levels.size());
  IntegerVector sample_sizes(sample_levels.size());
  std::vector<int> global_to_local(classes, -1);
  for (std::size_t sample_index = 0; sample_index < sample_levels.size(); ++sample_index) {
    const int sample = sample_levels[sample_index];
    checkUserInterrupt();
    const auto& rows = groups[sample - 1];
    sample_sizes[sample_index] = static_cast<int>(rows.size());
    std::vector<int> local_levels;
    local_levels.reserve(rows.size());
    for (int row : rows) local_levels.push_back(labels_r[row] - 1);
    std::sort(local_levels.begin(), local_levels.end());
    local_levels.erase(
      std::unique(local_levels.begin(), local_levels.end()), local_levels.end());
    for (std::size_t local = 0; local < local_levels.size(); ++local) {
      global_to_local[local_levels[local]] = static_cast<int>(local);
    }
    std::vector<int> local_labels(rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i) {
      local_labels[i] = global_to_local[labels_r[rows[i]] - 1];
    }
    const fiber_margin::CoordinateBlock coordinates =
      fiber_margin::robust_unit_coordinates(xy, rows);
    dimensions_used[sample_index] = coordinates.dimensions;
    const fiber_margin::SampleResult result = fiber_margin::refine_sample(
      coordinates.values, local_labels, coordinates.dimensions,
      static_cast<int>(local_levels.size()), settings);
    for (std::size_t i = 0; i < rows.size(); ++i) {
      const int row = rows[i];
      output[row] = local_levels[result.labels[i]] + 1;
      candidate[row] = local_levels[result.candidate[i]] + 1;
      margin_score[row] = result.margin_score[i];
      required[row] = result.required[i];
      dispersion[row] = result.dispersion[i];
      isolation[row] = result.isolation[i];
    }
  }
  return List::create(
    _["labels"] = output,
    _["candidate"] = candidate,
    _["margin_score"] = margin_score,
    _["required"] = required,
    _["atlas_dispersion"] = dispersion,
    _["isolation"] = isolation,
    _["changed"] = output != labels_r,
    _["sample_codes"] = wrap(sample_levels),
    _["dimensions_used"] = dimensions_used,
    _["sample_sizes"] = sample_sizes
  );
  END_RCPP
}
