---
title: "FiberMargin: Two-Sided Path Enclosures for Multiclass Spatial Label Repair"
author: "Moussa Kassim; Martin Ocharo; Dalia Ahmed; Dupe Ojo; Alessia Vignoli; Leonardo Tenori; Dinesh Gupta; Silvano Piazza; Stefano Cacciatore"
bibliography: jmlr/fibermargin.bib
link-citations: true
---

```{=latex}
\begin{abstract}
```

Multiclass spatial label fields can contain isolated errors, unstable borders, and
mislabeled patches. We study
training-free repair from coordinates and one imperfect categorical field, without
images, logits, anchors, or reference labels. FiberMargin is a deterministic
two-sided enclosure operator. Rotated Hilbert paths support forward and reverse
exponential recurrences; adding their square-root amplitudes rewards rival-class
evidence that encloses a query from both directions. A label changes only when the
resulting atlas contrast clears a barrier determined by path disagreement, class
prevalence, and local isolation. We derive an exact repair--damage decomposition,
characterize the enclosure field, and establish finite-input stability.

Against eight equal-input controls, FiberMargin obtained
the highest Accuracy and adjusted Rand index on 708 complete simulations, 108 unseen
geometries, procedural multiclass masks, DLPFC, and MERFISH, and the highest Accuracy
on variable-density 3D fields. It narrowly trailed multiscale mode on CRC. On saved
CamVid and S3DIS predictions it improved the upstream labels but did not lead the
comparison, delimiting its transfer range. The C++ implementation processed 500,000
observations in 6.15 seconds. The evidence supports FiberMargin as an auditable
operator for irregular multiclass spatial fields, not as a universal mask smoother.

```{=latex}
\end{abstract}
\begin{keywords}
multiclass spatial labels, categorical field repair, noisy labels, space-filling
curves, selective prediction
\end{keywords}
```

# Introduction

Semantic masks, anatomical label maps, spatial-omics domains, land-cover maps, and
point-cloud annotations are all categorical fields on a spatial domain. Once such a
field exists, a practitioner may trust its large-scale organization but suspect local
assignment errors. Re-running the upstream annotator or segmentation model may be
impossible, expensive, private, or unnecessarily disruptive.

We study *training-free multiclass spatial label repair*. The input is one label field
$\tilde y_i\in\{1,\ldots,K\}$ at coordinates $x_i\in\mathbb R^d$; the output is a
repaired field and a pointwise audit map. Source images, expression values, logits,
clean masks, and a training collection are unavailable. Regular pixel and voxel
arrays are supported, while irregular 2D/3D point sets are the primary setting. This
information budget separates repair from segmentation and from learning robustly
with noisy annotations. The package also contains a separate binary utility, evaluated
only as an auxiliary diagnostic in the supplement; it is not evidence for the
multiclass enclosure mechanism developed here.

Section 2 situates this coordinate-and-label-only setting among classical spatial
restoration, graphical models, learned segmentation refiners, noisy-label learning,
annotation fusion, self-supervised denoising, and topology-aware segmentation. The
supplement summarizes the information available to each method family.

FiberMargin represents each specimen through multiple one-dimensional charts made by
rotating the coordinates and ordering them with a Hilbert space-filling curve
[@sagan1994; @moon2001]. Each chart supports a linear-time recurrence. A
single fixed transport range propagates unit class evidence in both path directions.
The two directional amplitudes form an enclosure field that rewards bilateral support,
while large coordinate jumps attenuate transport. The central decision rule is simple:
accept a rival only when its mean rival-versus-observed contrast across the path atlas
exceeds the barrier induced by chart disagreement, class prevalence, and local query
isolation. The term *margin* refers to this difference between rival and observed-class
evidence.

The paper makes four empirical and theoretical contributions:

1. an information-budget taxonomy for categorical-field cleaning, a formal
   training-free repair interface, and an exact identity separating repaired
   errors from damage to correct labels;
2. a two-sided enclosure margin based on an atlas of chain orderings, with an exact
   amplitude characterization, conditional chart-interior recovery results, and
   finite-input stability;
3. a reproducible benchmark spanning procedural multiclass masks, spatial
   simulations, volumetric fields, controlled biological coordinates, and two real
   upstream-prediction evaluations, with calibrated Potts,
   conditional-random-field (CRF), propagation, and multiscale-mode controls under
   the same information budget;
4. a minimal R API and C++ implementation measured through 500,000 observations and
   20 classes.

# Related Work and Task Positioning

## Classical categorical-field restoration

Cleaning a spatial label field predates modern segmentation networks. Mathematical
morphology provides erosion, dilation, opening, closing, reconstruction, and
connected filtering driven by a chosen structuring element [@serra1982; @soille2003].
Mode filters and connected-component rules similarly encode a local scale or an
explicit size threshold. These operators are transparent and training free, but their
behavior is tied to a lattice neighborhood and a fixed geometric template; applying
the same template to a rare island, a narrow class, and an isolated error can produce
opposite outcomes.

Conditional random fields (CRFs) are graphical models that combine a unary preference
for each class with pairwise compatibility terms between locations [@geman1984;
@besag1986]. Iterated conditional modes is a fast local optimizer, while graph-cut
constructions give strong optimization results for important discrete energies
[@boykov2001]. Fully connected CRFs add long-range pairwise interactions and typically
combine positions with color or learned unary scores [@krahenbuhl2011]. CRF-RNN
(*conditional random fields as recurrent neural networks*) unrolls CRF mean-field
updates as recurrent layers, allowing them to be trained with a segmentation network
[@zheng2015crfrnn]. Random-walker segmentation and harmonic label propagation spread
information across an affinity graph [@grady2006; @zhu2003]. FiberMargin is a
coordinate-and-label-only restoration operator: it treats every observed label as
potentially erroneous, rather than designating trusted seeds or using calibrated class
probabilities. It aggregates spatial evidence along several rotated Hilbert paths, each
evaluated by a one-dimensional recurrence. Under this same input constraint, we compare
fixed Potts ICM and alpha-expansion Potts inference, a coordinate-only CRF,
self-anchored harmonic propagation, and multiscale mode filtering. Methods that also
use image features, calibrated unaries, or trusted seeds operate in richer input
settings and are discussed separately.

## Feature-guided and learned mask refinement

A large computer-vision literature improves an existing segmentation by returning to
the source image or its learned representation. GrabCut couples user initialization,
color models, and graph-cut optimization [@rother2004]. PointRend adaptively resamples
uncertain points using backbone features [@kirillov2020pointrend]. CascadePSP refines
a coarse mask from global and local image context [@cheng2020cascadepsp], SegFix learns
directions from boundary pixels toward reliable interiors [@yuan2020segfix], and
Boundary Patch Refinement retrains high-resolution image-and-mask patches around
predicted contours [@tang2021bpr]. Mask Transfiner detects incoherent quadtree regions
and updates them with learned image features [@ke2022masktransfiner]. These methods
principally target contour fidelity or high-resolution instance masks.

More recent systems learn a generic repair prior rather than attaching refinement to
one upstream architecture. SegRefiner learns a conditional discrete diffusion process
from coarse/clean mask pairs [@wang2023segrefiner], while RNCA learns recurrent
self-repair using image context and reference masks [@silbernagel2026rnca]. Both learn
from paired masks, and RNCA additionally uses image context. They address a
training-enabled setting that is distinct from the single-field repair task studied
here. Our external-mask experiments test the stricter mask-only operator, not
superiority over feature-guided refiners.

## Noisy-label learning, dataset cleaning, and blind denoising

Learning with noisy labels usually asks how to estimate a predictor from many
feature--label examples despite corruption. Representative strategies estimate a
noise transition and correct the loss [@patrini2017], exchange small-loss examples
between networks [@han2018coteaching], infer likely dataset label errors from
out-of-sample class probabilities [@northcutt2021], or use instance-level confidence
information [@berthon2021]. Their output is a robust predictor, a ranked error list, or
a revised training set. FiberMargin does not train a predictor. It revises one
supplied categorical field once, using only its coordinates and observed labels.

Pixelwise annotation noise has also motivated segmentation-specific training methods.
ADELE uses early network predictions to update noisy masks during training
[@liu2022adele]; Spatial Correction models correlated annotation errors and
progressively corrects them [@yao2023spatialcorrection]; and PNAL combines prediction
history with clusterwise correction for noisy 3D point-cloud training
[@ye2021pointnoise]. These methods are highly relevant to the application but require
images or point features, a trainable segmentation model, and repeated optimization.
They address *learning from noisy masks*, whereas our task is *repairing a supplied
mask when the learning system is unavailable*.

Noise2Self provides a closer conceptual analogy: leave-out prediction can denoise
measurements without clean targets when noise is conditionally independent across
measurement dimensions [@batson2019noise2self]. FiberMargin's leave-self-out transport
likewise prevents a query label from certifying itself. The objectives nevertheless
differ. Noise2Self estimates or trains a denoising function for numerical measurements
under a $\mathcal J$-invariance argument; FiberMargin makes categorical edits from
spatial redundancy and does not estimate deployment risk from the noisy field alone.

## Repeated annotations, topology, and selective action

When several labels of the same object are available, aggregation is statistically
better identified than one-mask repair. Dawid--Skene estimates annotator confusion
rates [@dawid1979], and STAPLE jointly estimates a latent segmentation and the
performance of multiple raters or algorithms [@warfield2004]. Multi-annotator
segmentation models extend this idea using images and learned annotator models
[@tanno2019annotators; @zhang2023annotators]. FiberMargin takes one label field as
input, so it cannot compare annotators or estimate how reliable each one is. When
several annotations of the same object are available, they should be combined with an
annotation-fusion method or studied in a separate multi-annotation analysis.

Topology-aware segmentation addresses a different failure mode. Persistent-homology
losses constrain Betti structure during network training [@hu2019topology], and clDice
targets connectivity of tubular structures through centerlines [@shit2021cldice].
These methods use clean topology-bearing references or an explicit structural target.
FiberMargin has neither. Its adaptive barrier is instead related to selective
prediction: it acts only when evidence clears a threshold [@elyaniv2010], but the
deterministic score is not a calibrated probability or a topology certificate. We
therefore report Boundary IoU [@cheng2021boundaryiou], rare-class measures, correction
recall, and damage, and retain the thin-network failure as a counterexample.

## Spatial omics and the exact comparison scope

Spatial transcriptomics is the principal irregular-domain application here. SpaGCN,
GraphST, BayesSpace, BANKSY, and KODAMA estimate spatial domains from gene expression
and spatial coordinates, with histology incorporated by some methods
[@hu2021; @long2023; @zhao2021; @singhal2024; @abdelshafy2025kodama]. FiberMargin starts
only after a categorical assignment exists and never claims to rediscover domains. We
compare the optional post-hoc correction rules from SpaGCN and GraphST, not their
clustering models, because those rules accept labels and coordinates. The same distinction
applies beyond transcriptomics: a point-cloud network trained under noisy labels and
a post-hoc repair of its exported semantic field are complementary stages.

Across these literatures, comparability is governed by six questions: Is there one
mask or repeated annotations? Are source measurements, image features, or logits
available? Are any labels trusted? Is a training collection available? Is topology
specified? Must the method work natively on irregular 2D/3D coordinates? FiberMargin
occupies the restrictive cell with one fallible field, coordinates only, no fitted
corpus, no anchors, and no topology template. The contribution is a new decision
operator for this cell, not a claim to have invented Hilbert orderings, exponential
attenuation, relative margins, or selective decisions [@rockafellar2000; @elyaniv2010].
Its defining multiclass construction is the two-sided path enclosure: opposite
directional states are combined before class normalization, and a rival must remain
positive across distinct spatial flattenings strongly enough to clear their observed
disagreement. The Supplement compares this operator with additive transport, one
reference path, smaller atlases, and reduced admission barriers.

# FiberMargin

## Mask-repair objective and specimen isolation

Let $y_i^*$ denote an evaluation-only reference, $\tilde y_i$ the imperfect input,
and $\hat y_i$ the repaired output. Define initial and repaired accuracy as $A_0$ and
$A_1$, correction recall $R=\Pr(\hat y=y^*\mid\tilde y\ne y^*)$, and damage
$D=\Pr(\hat y\ne y^*\mid\tilde y=y^*)$. Partitioning cells by initial correctness
gives the exact identity

$$
A_1=A_0(1-D)+(1-A_0)R,
\qquad
A_1-A_0=(1-A_0)R-A_0D.
$$

Thus mask cleaning is not synonymous with smoothness: an edit policy improves
accuracy only when corrected error mass exceeds damaged correct mass. We report both
terms, together with class-balanced and boundary measures, rather than selecting a
method by accuracy alone. At deployment $y^*$ is absent; FiberMargin therefore
returns a candidate, its uncalibrated `margin_score`, the required barrier,
`repair_margin = margin_score - required`, and the final retain/change decision.
These diagnostics are deterministic evidence contrasts, not probabilities or
estimates of expected loss.

Let $s_i\in\{1,\ldots,S\}$ identify a specimen or independent mask. FiberMargin applies the complete
operator separately to each set $I_s=\{i:s_i=s\}$. Coordinates, the locally observed
class alphabet, class counts, path states, and thresholds never cross specimens.
This is stronger than assigning zero weight to distant sections: the output and every
pointwise diagnostic for one specimen are mathematically invariant to all coordinates
and labels in every other specimen, including specimens that use the same numeric
coordinate range.

Within specimen $s$, axes with zero empirical range are removed before normalization.
Each remaining coordinate dimension is mapped by its empirical 1st and 99th percentiles
to a robust unit interval and clipped to $[-0.1,1.1]$. We write the normalized coordinate
as $z_i\in\mathbb R^{d_s}$, where $d_s\in\{2,3\}$; inputs with fewer than two varying
axes in any specimen are rejected. Thus a constant third coordinate invokes exactly the
2D operator, whereas a varying third coordinate participates in rotations, Hilbert
ordering, path distances, and the binary distance calculation.

## Path atlas and two-sided enclosure

For view $a$, a fixed rotation $R_a$ is applied to $z_i$. Each rotated coordinate
axis is rescaled to $[0,1]$ for a 16-bit Hilbert code, yielding a permutation
$\pi_a=(\pi_{a1},\ldots,\pi_{an})$. The atlas contains nine paths in 2D and 12 in
3D. Each path is determined solely by its rotation and Hilbert ordering.

Let $\Delta_{at}=\|z_{\pi_{at}}-z_{\pi_{a,t-1}}\|_2$ and let $r_a$ be the median
positive path step. At an interior query $i=\pi_a(t)$ define
$q_{ai}=\sqrt{\Delta_{a,t-1}\Delta_{at}}$; at an endpoint use the available step.
The transport action and query-isolation factor are

$$
h_{at}=\frac{\Delta_{at}}{r_a},\qquad
\omega_{ai}=\max\left\{1,\frac{q_{ai}}{r_a}\right\}.
$$

For class indicator $z_{ajc}=\mathbf 1\{\tilde y_j=c\}$, one forward and one reverse
recurrence at the fixed range $\lambda=5$ give

$$
L_{aic}=\sum_{j\prec_a i}e^{-D_a^+(i,j)/5}z_{ajc},\qquad
R_{aic}=\sum_{j\succ_a i}e^{-D_a^-(i,j)/5}z_{ajc},
$$

where $D_a^+$ and $D_a^-$ are sums of the directed actions. The strict path
inequalities are leave-self-out: a query's supplied label is inserted only after its
directional state has been read. FiberMargin combines directional amplitudes as

$$
H_{aic}=(\sqrt{L_{aic}}+\sqrt{R_{aic}})^2,\qquad
p_{aic}=\frac{H_{aic}}{\max\{\sum_r H_{air},10^{-12}\}}.
$$

The cross term rewards evidence that encloses a query from both path directions while
retaining unilateral evidence. Sources have unit mass; local isolation affects only
the final acceptance decision.

## Atlas margin and selective decision

Relative to the supplied class $\tilde y_i$, define the path contrast
$m_{aic}=p_{aic}-p_{ai,\tilde y_i}$. For $A$ paths, the atlas margin and descriptive
chart-disagreement scale are

$$
m_{ic}=\frac1A\sum_a m_{aic},\qquad
s_{ic}=\frac1A\left\{\sum_a(m_{aic}-m_{ic})^2\right\}^{1/2}.
$$

The second expression is the empirical population standard deviation across the
deterministic paths divided by $\sqrt A$; it is not a calibrated confidence interval.
The candidate $c_i^*=\arg\max_c m_{ic}$ is initialized at the supplied class, so exact
ties retain that class. With observed class count $n_c$, median positive class count
$n_{\mathrm{med}}$, and median isolation $\omega_i=\operatorname{median}_a\omega_{ai}$,

$$
\delta_c=\sqrt{\frac{n_{\mathrm{med}}}{\max(n_c,1)}},\qquad
B_{ic}=\frac{s_{ic}\omega_i}{\delta_c}.
$$

FiberMargin changes $\tilde y_i$ to $c_i^*$ only when $c_i^*\ne\tilde y_i$ and
$m_{i c_i^*}\ge B_{i c_i^*}$. The returned `margin_score` is $m_{i c_i^*}$,
`required` is $B_{i c_i^*}$, and `repair_margin` is their difference. These are
uncalibrated deterministic contrasts, not probabilities.

For exactly two observed classes, the implementation uses one strict local ballot
over $\min(n-1,2^{d+4})$ nearest other locations; a tie retains the supplied class
and a distinct candidate requires one net vote. Specimens with fewer than 20
locations, fewer than two observed classes, or no coordinate extent are returned
unchanged. In a binary field the rival is unique, so the signed difference between
local rival and supplied-class counts provides the complete local candidate
comparison without multiclass path normalization. This binary specialization is
fixed and is not user-selectable.

## Characterization and conditional guarantees

The square-root enclosure has a direct characterization. Within fields of the form
$F(L,R)=(g(L)+g(R))^2$ with nonnegative $g$ and $g(0)=0$, exact retention of
one-sided evidence, $F(L,0)=L$ and $F(0,R)=R$, forces $g(u)=\sqrt u$ and therefore
$F(L,R)=(\sqrt L+\sqrt R)^2$. This applies to the unnormalized field and does not
derive the atlas size, rotations, or range.

The Supplement proves deterministic path-weight, finite-input recurrence-stability,
and specimen-isolation results for the implemented recurrences. Its binary and
multiclass recovery theorems are conditional chart-interior results under stated
witness and label-noise assumptions, not universal consistency claims. In particular,
they do not imply recovery of coherent regional overwrites or correctness at every
boundary.

# Experimental Design

The evaluation separates six questions that a single headline accuracy cannot
answer: whether the operator corrects known errors under controlled geometry and
density; whether it transfers to unseen shapes and native 3D coordinates; whether it
works on real tissue arrangements without expression or images; whether a regular
raster instead favors a grid-specific smoother; which parts of the operator matter;
and whether the C++ implementation remains practical at large sample sizes. In every
corruption study, the reference field is retained only for scoring. It is never an
input to FiberMargin or a comparator.

## Evaluation units

For simulation inference, the independent unit is a geometry family rather than one
corruption draw: repeats on the same field share its class layout, boundary geometry,
and sampling pattern. Confidence intervals therefore use 10,000 paired bootstrap
resamples of whole-family means. Real-coordinate perturbations are descriptive
recovery tests on a fixed tissue arrangement, not biological replicates.

For the genuine upstream-mask evaluation, an independently annotated CamVid test frame
is the summary unit. Frames retain sequence dependence, so image-wise summaries are
descriptive rather than formal out-of-domain significance tests.

## Controlled spatial recovery studies

**Purpose.** Controlled fields make the correct class known at every location, so they
can measure both repaired errors and newly damaged correct labels while varying one
spatial difficulty at a time. The complete simulation suite contains 708 cases from 12
planar or layered geometry families, including jagged and curved stripes, thin layers,
branching regions, rings, disconnected components, density gradients, 2D/3D A-B-C
gradients, five region-concentration profiles, four corruption mechanisms, and one or
three specimens. It asks whether recovery persists across the forms of class shape,
sampling imbalance, and local error that occur in spatial annotations. A separate
108-case unseen-shape suite uses six held-out geometry families; it tests transfer to
new formulas rather than another draw from the primary shapes. Every procedural
corruption retains one correct exemplar of each true class represented in a specimen.
The target is repair of labels that remain represented in the supplied field, rather
than discovery of a class absent from that field.
Figure 1 compares Accuracy and ARI across both simulation ledgers and the dedicated
variable-density 3D panel.

The separate, non-overlapping volumetric suite asks whether the same operator works
on genuinely three-dimensional arrangements rather than stacked planar examples. It crosses eight
3D structures with uniform, class-imbalanced, or irregular-plane acquisition and
random, boundary, patch, or regional corruption. It contains 288 cases of 12,000
points and 32 size-confirmation cases of 60,000 points.

Controlled biological-coordinate studies ask the same question on real spatial
arrangements while keeping recovery measurable: annotations are first treated as the
reference and then locally corrupted. This design does not claim correction of
unadjudicated biological annotation errors. We use 47,329 spots from 12 DLPFC sections
and three donors [@maynard2021; @pardo2022] to test layered, multi-section tissue;
28,317 cells from five consecutive MERFISH hypothalamic planes from one mouse
[@moffitt2018; @li2022bass] to test a registered multi-plane arrangement; and 194,541
annotated locations from a 19-class colorectal cancer (CRC) Visium HD tissue derived
from the 10x Genomics *Visium HD Spatial Gene Expression Library, Human Colorectal
Cancer (FFPE)* dataset, analyzed with Space Ranger 4.0.1 [@tenxcrc2025], to test a
dense, multiclass tissue field. The CRC WSI region classes were derived by the
authors. Labels are replaced by labels from plausible adjacent regions under fixed
random, interface, compact-patch, and regional protocols. FiberMargin receives only
the resulting coordinates, corrupted labels, and specimen identifiers.

## Procedural raster-mask diagnostic

This diagnostic asks a deliberately difficult question for a general coordinate-based
operator: can it repair a regular categorical raster when no source image or model is
available, and how much does a fixed grid-local smoothing prior help? The suite crosses
six 96x96 families--nested regions, jagged layers, thin networks, disconnected
islands, touching regions, and perforated rings--with three geometry seeds,
independent, boundary, or coherent-patch errors, and 15% or 30% corruption, for 108
cases. The families respectively stress containment, irregular boundaries, narrow
connectivity, rare components, close class contacts, and holes. Every method receives
only the corrupted categorical array. The same exported package functions create
corruptions, clean masks, and score outputs.

## Auxiliary binary building-mask diagnostic

This is a controlled cross-domain geometry test. It asks whether the unchanged
mask-only operator can repair held-out non-transcriptomic mask geometry, while
controlled corruptions make the exact recovery and damage observable. We used 55 256x256 masks
from a held-out project-separated test split of the HOT Very-High-Resolution Building
Segmentation Dataset [@hot2026buildings]. Eligibility required at least 256 pixels of
each class. Only binary masks were read, never the aerial images. Each source mask
received random, boundary-concentrated, and compact-patch corruption at 10% and 25%,
producing 330 cases. This tests recovery from injected defects, not segmentation from
imagery or correction of naturally occurring mapping errors.

Because every HOT field is binary, this study evaluates the package's fixed local
ballot utility rather than the multiclass path-enclosure operator. Its results are
reported only in the supplement as a package-scope diagnostic and are not used as
evidence for the enclosure mechanism.

The matched comparison contains nine coordinate-and-label rules: FiberMargin, fixed
Potts ICM, alpha-expansion optimization of a fixed Potts energy, a coordinate-only
CRF, self-anchored harmonic propagation, multiscale mode filtering, fixed SpaGCN
correction and GraphST correction rules, and a one-pass 8-NN modal filter. Every rule receives the
same coordinates, observed categorical field, and specimen identifier in every study.
Results average the six corruption conditions within each source mask, which is the
independent evaluation unit.

## Controlled public-annotation diagnostics

To test ordinary segmentation-mask geometry without turning the study into an
image-refinement benchmark, we used the categorical annotations of 24 Oxford-IIIT Pet
images [@parkhi2012] and 24 CamVid images [@brostow2009]. Oxford trimaps were converted
to foreground/background and CamVid retained its annotated semantic classes, with void
pixels excluded. Each mask was nearest-neighbor resized to a maximum side length of
128, then subjected to random, boundary, and compact-patch errors at 15% and 30%, for
144 corruptions per collection. RGB imagery was never read. These are controlled
annotation-corruption diagnostics, not evaluations of naturally predicted masks.

## Real upstream-mask evaluations

We separately evaluate repair of unmodified upstream categorical predictions against
independent references. For CamVid [@brostow2009], a fixed Extremely Randomized Trees
semantic pixel classifier [@geurts2006extratrees] was fitted only on the official
training frames using RGB, locally averaged RGB, edge magnitude, and normalized pixel
position. It then predicted all 233 official test frames; the 11 retained semantic
classes were evaluated at a maximum side length of 160 with void reference pixels
excluded. The test annotations were never used to fit the upstream classifier.

For S3DIS [@armeni2016], the same classifier family was fitted to Areas 1, 2, 3, 4,
and 6 using XYZ and RGB, then evaluated on one deterministic 4,096-point block from
each of 68 Area 5 rooms. FiberMargin and the controls received only XYZ and the saved
categorical prediction. These deliberately simple upstream models create independent
prediction/reference pairs for testing the repair interface; they are not intended as
state-of-the-art segmentation benchmarks.

Images and point features are used only to produce the upstream fields. FiberMargin
and every equal-input control receive only a categorical prediction and coordinates,
never RGB values, detector scores, logits, or reference labels.

## Fixed comparator calibration

The four calibrated controls were set before the reported evaluations. A
non-overlapping development panel of 24 simulated 2D/3D fields used distinct geometry
formulae, corruption draws, and random seeds from the reported suites. For each
control family, we selected the setting with the largest mean accuracy subject to
the same absolute damage cap of 0.05; ARI broke accuracy ties. SpaGCN and GraphST retain their fixed published correction
settings, and the modal reference retains $k=8$. No comparator setting was selected
from a reported test result. The calibration files include every candidate setting and
case-level result.

## Comparators and metrics

The direct comparison asks whether a rule improves the same supplied categorical field
under the same information budget. Numerical eligibility therefore requires a post-hoc
categorical rule that consumes exactly one categorical field, coordinates, and specimen
identifiers, without images, expression, logits, trusted anchors, repeated annotations,
or a fitted training collection. Let $Y$ be the one-hot encoding of the observed labels
and let $P_k$ be a row-normalized $k$-nearest-neighbor coordinate matrix constructed
separately in each specimen. The four calibrated control families use only these
inputs.

The first family is a fixed local Potts-style ICM rule, and the second is standard
multi-label alpha-expansion graph-cut inference [@boykov2001] for

$$
E(z) = 5\sum_i \mathbf{1}\{z_i\ne\widetilde y_i\} +
\sum_{\{i,j\}\in G_8}\mathbf{1}\{z_i\ne z_j\},
$$

where $G_8$ is the undirected, symmetrized eight-nearest-neighbor graph built within
each specimen or mask. Alpha expansion starts from the observed field and performs
three complete class cycles. Its graph construction, unary weight, pairwise weight,
and cycle count are fixed across evaluations; the ICM rule is a separately configured
eight-pass local optimizer. These controls distinguish global discrete-energy
optimization from local smoothness updates.

The coordinate-only CRF is a sparse Gaussian-coordinate conditional random field
evaluated by eight mean-field updates,

$$
Q_i^{(t+1)}(c)\propto
\exp\left[-1.5\mathbf1\{c\ne\widetilde y_i\}
+2.5(P_{16}Q^{(t)})_{ic}\right].
$$

It uses no image features, logits, learned unaries, or reference labels. The
self-anchored harmonic-propagation control supplies every observed label as an equally
weighted soft observation, rather than designating a trusted seed set:

$$
F^{(t+1)}=0.20Y+0.80P_{24}F^{(t)},\qquad t=0,\ldots,11.
$$

Finally, the multiscale-mode control averages categorical support over
$k\in\{6,14,30\}$ coordinate neighborhoods, adds 0.20 support to the observed class,
and changes a label only when the best rival exceeds its observed-class support by
0.05. Thus it is a multiscale version of the mode-filter family, defined on the same
irregular 2D/3D coordinates as FiberMargin rather than only on a pixel lattice.
Alpha expansion, the coordinate-only CRF, self-anchored propagation, multiscale mode,
and Potts ICM are joined by the fixed SpaGCN correction rule, fixed GraphST correction
rule, and the one-pass 8-NN modal reference. These nine rules are evaluated on every
controlled and real-upstream benchmark: planar and volumetric simulations, CRC,
DLPFC, MERFISH, procedural masks, held-out building masks, and image masks. We compare
the coordinate-label correction rules rather than end-to-end clustering systems. The
two published-rule reproductions were checked pointwise against official source commits
[@spagcnsource; @graphstsource]. They matched all 9,600 predictions on continuous
coordinates; 5 of 648 lattice predictions differed because exact-distance ties have no
specified stable order.

The modal reference is a one-pass instance of classical local mode filtering
[@serra1982; @soille2003]. For each labelled location $i$, let
$\mathcal N_k(i)$ be its $k$ nearest non-focal labelled locations within the same
specimen or mask, and let $\widetilde y$ denote the corrupted labels. The filter returns
$$
M_k(i)\in\arg\max_{c}\sum_{j\in\mathcal N_k(i)}
\mathbf 1\{\widetilde y_j=c\}.
$$

The implementation resolves a tied mode with the nearest tied neighbour. The compact
$k=8$ rule is part of the matched nine-method roster. A broad $k=48$ version is a
separate regular-grid diagnostic, shown only for raster-mask studies; neither value is
tuned per mask or corruption condition. ADELE and PNAL require a trainable predictor,
data features, and prediction history; Noise2Self requires a numerical
self-supervision setting; image-guided refiners require image features or logits; and
annotator aggregation requires repeated labels. These methods are relevant neighboring
tasks, but their information contracts differ. Supplement S1 gives the full input
taxonomy.

No single score captures useful repair. The evaluation maintains four core measures:
accuracy and ARI measure global recovery, correction recall measures recovery among
initially wrong labels, and damage measures mistakes introduced among initially correct
labels. Each main-text table shows the subset most relevant to its question, while the
supplement reports the complete core outcomes and secondary diagnostic
measures: mean IoU for mask overlap; one-grid-step Boundary IoU
[@cheng2021boundaryiou] for geometric edge recovery; macro and worst recall for class
balance; and boundary and rare/sparse-region accuracy for difficult locations. For a
reference field $y^*$ and evaluation neighborhood $\mathcal N_i$, the generic boundary indicator is
$b_i=\mathbf 1\{\exists j\in\mathcal N_i:y_j^*\ne y_i^*\}$, with anatomical
class-pair restrictions for DLPFC and MERFISH. If $U(i)$ is an evaluation unit and
$r_i$ a reference region or class, let
$n_{Uc}=\sum_{i:U(i)=U}\mathbf 1\{r_i=c\}$,
$c_U^-=\arg\min_{c:n_{Uc}>0}n_{Uc}$, and
$s_i=\mathbf 1\{r_i=c_{U(i)}^-\}$; no prevalence-percentage cutoff is used. The
neighborhood sizes, the 20% geometric-boundary threshold used by simulations,
tie handling, and all dataset-specific exceptions are fixed in Supplement S4.
Correction recall and damage are the two conditional terms in the exact identity
above. The shared core metrics make the correction--protection trade-off visible in
every main-text comparison.

The scaling study is an engineering test of the complete R-to-C++ call rather than a
claim of a universal runtime ranking. Supplement S3 reports matched one-worker timing.
A timed call includes validation, specimen splitting, coordinate normalization,
geometry construction, inference, and output assembly; it excludes file I/O,
corruption generation, reference-derived masks, scoring, and plotting. Hardware,
compiler flags, and separate deterministic CPU-parallel verification are reported
there.

# Results

## Broad multiclass controlled recovery

**Question.** Does the two-sided enclosure margin improve multiclass labels across
prespecified and previously unseen spatial structures? The complete suite varies
geometry, density, class balance, dimensionality, and corruption mechanism; the
unseen suite introduces new geometry formulas; and the 3D panel tests native
variable-density volumes. Figure 1 summarizes 36 variable-density 3D fields, 708
complete simulations, and 108 unseen geometries. Every method receives the same
coordinates and corrupted categorical labels.

FiberMargin obtained the highest Accuracy and ARI on both the complete and
unseen-geometry ledgers. It also obtained the highest Accuracy on variable-density
3D, while multiscale mode had a 0.0010 higher ARI. Descriptive paired
family-clustered intervals compare FiberMargin with the strongest displayed control;
they are not selection-adjusted. The complete-suite Accuracy difference was 0.0112
(95% interval 0.0076--0.0146), and the unseen-suite difference was 0.0108
(0.0078--0.0147). In the variable-density 3D panel, the Accuracy interval crossed
zero and the ARI difference was compatible with a tie. Thus the strongest evidence
supports broad multiclass recovery across irregular geometries, not a universal
ranking over every spatial repair problem.

![Accuracy and adjusted Rand index across 36 variable-density 3D fields, 708 complete simulations, and 108 unseen geometries. These are evaluations of FiberMargin's multiclass two-sided enclosure mechanism. Every method receives the same coordinates and corrupted categorical labels; color and shape jointly identify the method.](figures/fibermargin_simulation_scope.png){width=100%}

## Controlled biological-coordinate recovery

This section reports three real-coordinate studies. **Question.** Given the actual
coordinate layout of layered tissue, consecutive planes, or a dense 19-class field,
can a method repair deliberately introduced local label errors without using
expression data, images, or the reference annotation? The published or author-derived
labels are treated as evaluation-only references, fixed neighbor-plausible
corruptions are introduced, and the original labels are hidden before inference.
These studies therefore establish controlled recovery on biological geometry, not
correction of naturally disputed annotations.

FiberMargin had the highest aggregate Accuracy and ARI on DLPFC
(0.8942 / 0.8220) and MERFISH (0.8905 / 0.8004) among the nine matched refiners.
CRC was the exception: multiscale mode narrowly led FiberMargin
(0.8604 / 0.7827 versus 0.8589 / 0.7806). FiberMargin recovered a larger fraction
of corrupted CRC labels (0.4444 versus 0.4105) but incurred more damage
(0.0478 versus 0.0336), making this a clear operating-point trade-off rather than a
uniform win. Detailed DLPFC and MERFISH operating points, boundary Accuracy,
sparse-region Accuracy, and class-balance comparisons are reported in the supplement.

The predefined random 25% CRC visualization illustrates what aggregate means can and
cannot show (Figure 2). FiberMargin raised Accuracy from 0.7500 to 0.9715 and
corrected 0.9204 of injected errors with damage 0.0114. In contrast, one 25%
coherent-patch case fell from 0.7500 to 0.7104. A contiguous patch overwritten with
a neighboring class is observationally ambiguous when only coordinates and the
categorical field are available; the failure is retained in the scenario archive.

![Full-tissue and zoomed 19-class CRC Visium HD recovery under 25% random corruption. FiberMargin uses its multiclass two-sided enclosure mechanism. The black rectangle identifies the zoom.](../benchmarks/results/fibermargin_final/colorectal_19class_full_and_zoom_indexed.png){width=100%}

## Procedural multiclass mask repair

**Question.** Does the same multiclass enclosure remain useful on regular categorical
rasters, where grid-specific local filters have a natural advantage? The suite crosses
six 96-by-96 geometry families, three error mechanisms, two corruption levels, and
three seeds, producing 108 cases. Table 1 includes the complete fixed comparator
roster plus a separate 48-NN regular-grid diagnostic.

The 48-NN filter had the highest control Accuracy (0.9070), while the 8-NN filter
had the highest control ARI and correction recall. FiberMargin reached 0.9113
Accuracy and 0.7392 ARI, leading both primary scores, and corrected 0.6210 of altered
locations. Its 0.0056 damage was higher than multiscale mode (0.0011) and alpha
expansion (0.0009), so the result is an Accuracy--damage trade-off rather than a free
gain from changing more cells.

| Method | Accuracy | ARI | Correction | Damage |
|:--|--:|--:|--:|--:|
| Initial | 0.7750 | 0.4242 | 0 | 0 |
| SpaGCN correction | 0.8740 | 0.6335 | 0.4573 | 0.0018 |
| GraphST correction | 0.9059 | 0.6964 | 0.6362 | 0.0163 |
| 8-NN modal filter | 0.9023 | 0.7218 | **0.6456** | 0.0231 |
| 48-NN modal filter | 0.9070 | 0.7012 | 0.6379 | 0.0153 |
| Alpha-expansion Potts | 0.8997 | 0.7037 | 0.5511 | 0.0009 |
| Coordinate CRF | 0.8662 | 0.6174 | 0.4166 | **0.0001** |
| Harmonic propagation | 0.8961 | 0.6913 | 0.5357 | 0.0019 |
| Multiscale mode | 0.9043 | 0.7210 | 0.5781 | 0.0011 |
| Potts ICM | 0.8988 | 0.7009 | 0.5498 | 0.0017 |
| FiberMargin | **0.9113** | **0.7392** | 0.6210 | 0.0056 |

: Aggregate core metrics on 108 procedural multiclass mask-repair cases. Higher is
  better except for damage. Every matched comparator receives the same corrupted
  categorical field and coordinates; the 48-NN modal filter is an additional
  regular-grid diagnostic. Secondary overlap, boundary, and rare-class metrics are
  reported in the supplement.

The operating-point audit multiplied the fixed barrier by two, retained the default,
or removed the barrier. Correction recall and damage were respectively 0.584/0.0021,
0.621/0.0056, and 0.645/0.0149. Qualitatively, FiberMargin cleaned isolated and
disconnected regions but removed much of one very thin network. That counterexample,
retained in the supplement, prevents a topology-preservation claim.

## Operator diagnostics

These diagnostics test the final operator rather than add method comparisons. The
supplement reports the complete atlas-size audit with paired Accuracy, ARI,
worst-recall, sparse-region, and damage outcomes. The results support the fixed
9-path 2D and 12-path 3D atlas used by the package; changing CPU workers does not
change the mathematical operator.

The specimen-isolation test verifies the stated mathematical contract rather than
recovery Accuracy: translating, rotating, coordinate-permuting, or label-permuting a
distractor specimen changed zero of 2,600 target predictions. This exact invariance
follows from independent specimen processing.

The software contract additionally tests exact equality of all pointwise diagnostics
under one, two, and four requested CPU workers; row-permutation equivariance;
translation and uniform-scale invariance; equality of 2D and constant-$z$ inputs; and
non-equivalence after removing informative $z$ from a genuinely 3D field. A mixed call
may therefore contain both effective 2D and effective 3D specimens without sharing
coordinate normalization or class evidence.

## Scaling

This feasibility study measures the complete R-to-C++ call rather than a universal
runtime ranking. Fresh R processes measured complete 20-class calls at 10,000,
50,000, 100,000, 250,000, and 500,000 observations. Runtime rose from 0.067 to
6.146 seconds. The fitted log-log slope was 1.166 ($R^2=0.998$). At 500,000
observations, peak resident memory was 569 MB and the increase during refinement was
422 MB. For fixed atlas size, transport range, and class count, time is
$O(n\log n+nK)$ and stored chart actions require $O(nK|\mathcal A|)$ memory.
Matched cross-method timings are reported in Supplement S15.

## External transfer and auxiliary scope tests

The final evaluations ask whether the operator transfers beyond controlled recovery
to categorical masks produced by real upstream predictors. All repair methods receive
the same saved categorical prediction and coordinates; source features and upstream
outputs beyond the categorical prediction are withheld. Figure 3 shows the matched
Accuracy--damage operating points, and the supplement reports the complete image- and
room-level matrices.

FiberMargin improved CamVid Accuracy by 0.0186 (descriptive 95% image-bootstrap
interval 0.0175--0.0197) and S3DIS Accuracy by 0.0204 (room-block interval
0.0152--0.0259). GraphST correction had the highest Accuracy in both studies; its
advantage over FiberMargin was 0.0061 on CamVid and 0.0112 on S3DIS. Harmonic
propagation also provided a more favorable Accuracy--damage point in these two
settings. Frames and room blocks retain sequence or scene dependence, so the
intervals are descriptive. These tests establish that FiberMargin can improve saved
predictions, but they also delimit its current transfer range: the strongest evidence
remains irregular multiclass controlled recovery.

![Accuracy gain and correction damage for saved CamVid and S3DIS upstream predictions. Both evaluations are multiclass and use FiberMargin's two-sided enclosure mechanism. Every repair method receives only the categorical prediction and coordinates; color and shape jointly identify methods.](figures/fibermargin_external_natural_accuracy_damage.png){width=85%}

The supplement separately reports 330 controlled corruptions of 55 held-out binary
building masks. That study evaluates the package's local-ballot utility, not the
multiclass enclosure developed in this paper, and is therefore excluded from the
primary evidence hierarchy.

# Discussion

FiberMargin addresses a deliberately restricted multiclass repair task with one
explicit margin-generating recurrence. Each remote observation contributes through a
known exponential path weight, the fixed path atlas determines the contrast, and a
positive gap makes the decision stable. The C++ implementation retains this structure
closely enough for source-level audit.

The experiments identify a coherent regime of advantage. FiberMargin led the complete
and unseen-geometry simulation ledgers, procedural multiclass masks, DLPFC, and
MERFISH, and attained the highest Accuracy in the variable-density 3D study. These
settings share multiclass structure, irregular or heterogeneous spatial support, and
errors that remain locally contradicted from more than one direction. CRC exposed an
Accuracy--damage trade-off rather than a clear loss: multiscale mode narrowly led the
primary scores, while FiberMargin corrected more corrupted labels and damaged more
initially correct labels.

The external-prediction studies define the boundary of that evidence. FiberMargin
improved both saved upstream predictions, but harmonic propagation and GraphST
correction provided stronger operating points. The auxiliary binary ballot was also
not the strongest rule on held-out building masks. These outcomes do not negate the
multiclass enclosure result; they show that regular binary masks and errors inherited
from an upstream predictor can favor simpler local operators. Recovery and protection
must therefore be reported jointly, and method choice should follow the geometry and
error regime rather than a universal ranking.

## Scope of the mask-repair task

The same input contract covers regular image masks, 3D voxel labels, irregular point
clouds, and categorical maps, but the present evidence is strongest for multiclass
irregular-coordinate fields. It is most relevant when an upstream label field exists
but its image, logits, training pipeline, or domain-specific model cannot be used.
Potential tests include LiDAR labels [@ye2021pointnoise], land-cover maps
[@isprssemantic], and multi-annotator medical masks [@zhang2023annotators]. The
supplement prespecifies suitable datasets and guards.

The taxonomy in Section 2 also clarifies how future comparisons should be organized.
Single-mask, coordinate-only operators can be compared directly at equal information.
Feature-guided refiners should form a second track that receives the image or point
features; noisy-label methods should form a training-time track; and STAPLE-like
methods should form a repeated-annotation track. A topology-aware track should report
connectivity metrics and receive the structural information its methods require.
SegRefiner and RNCA can learn visual or anatomical repair priors from data;
FiberMargin cannot. Conversely, FiberMargin requires no fitted weights and applies
unchanged to irregular 2D/3D coordinates. Collapsing these settings into one ranking
would measure access to information as much as repair quality.

## Limitations and next steps

The procedural masks, auxiliary HOT study, and public-annotation diagnostics use
injected defects. CamVid and S3DIS instead supply genuine upstream predictions with independent
references, but two simple upstream models cannot establish uniform performance across
natural error mechanisms. The thin-network and coherent-patch failures show that
the operator is not topology preserving. A contiguous regional overwrite can remain
observationally ambiguous from geometry and one categorical field alone. The
conditional recovery theorems require chart-interior witnesses and independent
bounded label noise; they do not establish universal consistency at boundaries or
under structured regional corruption. The
finite atlas is not exactly rotation invariant.

The biological tests also remain controlled recovery studies. DLPFC contains several
donors, but CRC and MERFISH each contribute one biological unit; natural
annotation uncertainty was not independently adjudicated. Strong next tests require
independent tissues, additional paired imperfect/reference masks from upstream models,
and 3D point clouds. A feature-enabled veto should be evaluated separately from the
present mask-only operator. For current evidence, FiberMargin is a fast, auditable
multiclass spatial-label repair method whose strongest support is irregular
controlled recovery; external-prediction generalization remains an open problem.

# Software and Reproducibility

The R package exposes `clean_categorical_mask(mask, samples = NULL, workers = NULL)`
for pixel/voxel arrays and
`refine_spatial_labels(xy, labels, samples = NULL, workers = NULL)` for irregular
coordinates. `corrupt_categorical_mask()` and `evaluate_mask_cleaning()` implement
the public benchmark contract and exact repair--damage audit. The mathematical
operator, including path construction, directional recurrences, margin aggregation, and
decisions, is implemented in C++17. R validates and encodes inputs, allocates one total
CPU budget, and restores output types and row order. One native call performs
per-specimen axis reduction, coordinate normalization, local class encoding, and
refinement. Specimens reuse the same budget and native stages do not nest thread pools.
Independent charts and row ranges use deterministic standard C++ threads on macOS,
Linux, and Windows. The returned factor retains input levels and identifiers and carries
the seven pointwise diagnostics together with per-specimen dimensionality, sample size,
edit count, edit fraction, before/after class sets, and explicit removed-class summaries.
The implementation does not constrain the number of output classes. The repository contains the
versioned source, simulators, evaluators, manifests, scenario metrics, figures, and all
scripts used here. Compact DLPFC and CRC coordinate/annotation derivatives are
redistributed under their source terms. The CRC derivative is CC BY 4.0, records
attribution and modifications, and excludes expression counts and tissue imagery.
MERFISH source data are not redistributed; scripts document the official source and
license. The HOT and Oxford-IIIT Pet controlled studies use mask labels only. The
CamVid preparation script uses source images solely to create fixed upstream categorical
predictions; the R repair evaluator reads only those saved masks, coordinates, and
independent references, never image values or model scores. Its manifest records source
terms and the exact upstream recipe.

```{=latex}
\acks{The authors thank 10x Genomics and the maintainers of spatialLIBD, SpaGCN, GraphST, and the Humanitarian OpenStreetMap Team for making data or source code available for reproducible evaluation.\par\textbf{Funding.} The authors received no specific funding for this work.\par\textbf{Competing interests.} The authors declare that they have no competing interests.}
```
