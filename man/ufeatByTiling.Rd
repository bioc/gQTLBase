\name{ufeatByTiling}
\alias{ufeatByTiling}
\alias{balancedFeatList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
split featurenames of SummarizedExperiment according to tiling,
or to achieve simple balance within seqnames
}
\description{
split featurenames of SummarizedExperiment according to tiling,
without redundancies
}
\usage{
ufeatByTiling(se, tiling, maxlen=20)
balancedFeatList(se, maxlen=20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{se}{
instance of SummarizedExperiment
}
  \item{tiling}{
GRanges instance corresponding to a genomic tiling
}
  \item{maxlen}{
numeric ... list elements longer than maxlen are chopped up
to have this length, to foster load balancing
}
}
\details{
ufeatByTiling uses findOverlaps, balancedFeatList uses
split on seqnames and BBmisc::chunk
}
\value{
a list with elements of names(rowRanges(se)) corresponding
to the elements of the tiling
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
library(geuvPack)
data(geuFPKM)
library(Homo.sapiens)
au = paste0("chr", 1:22)
tg_500k = tileGenome(seqinfo(TxDb(Homo.sapiens))[au,], tilewidth=500000,
  cut.last.tile.in.chrom=TRUE)
sn = ufeatByTiling(geuFPKM, tg_500k)
summary(sapply(sn,length))
sn2 = balancedFeatList(geuFPKM)
summary(sapply(sn2,length))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
