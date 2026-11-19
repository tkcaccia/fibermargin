\name{evaluate_mask_cleaning}
\alias{evaluate_mask_cleaning}
\title{Evaluate categorical mask cleaning}
\usage{
evaluate_mask_cleaning(
  reference,
  initial,
  cleaned,
  elapsed = NA_real_,
  method = NULL
)
}
\arguments{
\item{reference}{Reference categorical matrix or three-dimensional array.}
\item{initial}{Imperfect input mask.}
\item{cleaned}{Output mask to evaluate.}
\item{elapsed}{Optional elapsed runtime in seconds.}
\item{method}{Optional method name included in the returned row.}
}
\value{A one-row data frame containing class-aware overlap, boundary agreement, rare-class, correction, and damage measures.}
\description{Evaluates a cleaned mask against a fixed reference and noisy input mask. The output includes overlap scores, class-boundary agreement, and correction/damage decomposition.}
\details{The evaluation domain consists of cells with non-missing labels in both the reference and the initial mask. A missing cleaned label on that domain is an error. Accuracy gain equals the repaired fraction minus the damaged fraction, both measured over the fixed domain.}
