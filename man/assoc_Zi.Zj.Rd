\name{assoc_Zi.Zj}
\alias{assoc_Zi.Zj}

\title{
Association between variable \eqn{Z_{i}} and \eqn{Z_{j}}.
}
\description{
Finding association between variables of \eqn{i^{th}}{i^th} position and \eqn{j^{th}}{j^th} position.
In any position wise aligned sequence dataset, occurences of R=(A,G) and Y=(C, T) at each position can be explained by a standard normal variate Z based on certain threshold value.
So, an association between any two position in the datast can be obtained which will be the association beween the two standard normal variate at this two positions.
However, the two nomal variates reprsenting the occurences of R and Y are independent of each other at a given position. 
}
\usage{
assoc_Zi.Zj(x)
}

\arguments{
  \item{x}{A dataframe of position wise aligned sequence dataset having A, T, G and C only.}
}
 
\details{
The user has to supply the sequence dataset in tab delimited format and not in FASTA format.
Each sequence (row) should contain only standard nucleotides (A, T, G and C).
Each sequence should be same length.
}

\value{
A numeric matrix of order \emph{L} by \emph{L} for the dataset of \emph{L} nucleotides long sequences.
}

\author{
Prabina Kumar Meher & A. R. Rao
}

\examples{
data(don_dat)
kk <- don_dat[1:300,]
zizj <- assoc_Zi.Zj(x=kk)
zizj
}
