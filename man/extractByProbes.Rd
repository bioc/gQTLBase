\name{extractByProbes}
\alias{extractByProbes}
\alias{extractByProbes,ciseStore,character,character-method}
\alias{extractByProbes,ciseStore,character,missing-method}
\alias{extractByRanges,ciseStore,GRanges-method}
\alias{extractBySymbols,ciseStore,character,character,character-method}
\alias{extractBySymbols,ciseStore,character,character,missing-method}
\alias{extractByRanges}
\alias{extractBySymbols}
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
extractBySymbols(store, symbols, sym2probe, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{store}{
instance of \code{\link{ciseStore-class}}
}
  \item{probeids, symbols}{
vector character tokens
}
  \item{gr}{
instance of \code{\link[GenomicRanges]{GRanges-class}}
}
  \item{sym2probe}{
named character vector of probeids with names given by corresponding
symbols
}
  \item{extractTag}{
character atom telling what field in the
archived GRanges is regarded as the probe or
gene identifier
}
  \item{\dots}{extra arguments to extractByProbes}
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
 store = makeGeuvStore()
 ebp = extractByProbes(store, c("ENSG00000183814.10", "ENSG00000174827.9"))
 ebp
 rr = range(ebp)
 ebr = extractByRanges(store, rr)
 ebr
 s2p = structure(c("ENSG00000183814.10", "ENSG00000163207.5", "ENSG00000228449.1", 
	"ENSG00000137962.8", "ENSG00000232848.1", "ENSG00000227280.1", 
	"ENSG00000238081.1", "ENSG00000117480.10", "ENSG00000253368.2", 
	"ENSG00000174827.9"), .Names = c("LIN9", "IVL", "RP11-177A2.4", 
	"ARHGAP29", "CTA-215D11.4", "RP11-458D21.2", "RP4-620F22.3", 
	"FAAH", "TRNP1", "PDZK1"))
 ss = extractBySymbols(store, c("IVL", "FAAH", "PDZK1"), s2p)
 ss
 }
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
