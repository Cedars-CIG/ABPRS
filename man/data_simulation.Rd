% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Simulation.R
\name{data_simulation}
\alias{data_simulation}
\title{Generate Simulation Data}
\usage{
data_simulation(m, effect_snps_vec, n, maf, effect_size, beta0, binary)
}
\arguments{
\item{m}{total number of snps}

\item{effect_snps_vec}{vector of length 4 containing the number of effective snps
with additive, dominant, recessive, and codominant encoding respectively}

\item{n}{total number of samples}

\item{maf}{minor allele frequency of all SNPs}

\item{effect_size}{effect size of the effective snps, which acts differently
depending on binary or continuous outcomes. See \eqn{\beta} in the equations sections.}

\item{beta0}{a value to ensure a certain percentage of case phenotype for
binary outcomes and to serve as a baseline measurement for continuous outcomes.
See \eqn{\beta_{0}} in the equations section.}

\item{binary}{a boolean flag indicating whether the data has binary outcomes (TRUE)
or continuous outcomes (FALSE).}
}
\value{
A dataframe containing the simulated phenotype in the first column and the
simulated genotype encoding in the rest of the columns.
}
\description{
This function generates simulated phenotype and genotype information.
}
\section{Equations}{


For binary outcomes, the phenotypes are simulated using a logistic model, as shown in the equation below: \cr
\deqn{\log \left(\frac{P(Y=1| \mathbf{X})}{1-P(Y=1| \mathbf{X})} \right)=\beta_0+\left(\mathbf{X}_{ADD}, \mathbf{X}_{Non-ADD} \right)\bm{\beta}}

For continuous outcomes, the phenotypes are simulated using the equation below: \cr
\deqn{{Y}=\beta_0+\left(\mathbf{X}_{ADD}, \mathbf{X}_{Non-ADD} \right)\bm{\beta}}
}

