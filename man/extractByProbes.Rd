\name{extractByProbes}
\alias{extractByProbes}
\alias{extractByRanges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
retrieve eqtlTest results from a ciseStore instance
}
\description{
retrieve eqtlTest results from a ciseStore instance
}
\usage{
extractByProbes(store, probeids, extractTag = "probeid")
extractByRanges(store, gr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{store}{
instance of \code{\link{ciseStore-class}}
}
  \item{probeids}{
vector character tokens
}
  \item{gr}{
instance of \code{\link[GenomicRanges]{GRanges-class}}
}
  \item{extractTag}{
character atom telling what field in the
archived GRanges is regarded as the probe or
gene identifier
}
}
\details{
an index will be searched if created by the ciseStore constructor
}
\value{
a GRanges instance
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%%%  ~~who you are~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
if (require(geuvStore)) {
 reg = partialRegistry()
 store = ciseStore(reg, addProbeMap=TRUE, addRangeMap=TRUE) 
 ebp = extractByProbes(store, c("ENSG00000183814.10", "ENSG00000174827.9"))
 ebp
 rr = range(ebp)
 ebr = extractByRanges(store, rr)
 ebr
 }
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
