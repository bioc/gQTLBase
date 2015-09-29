\name{describeStore}
\alias{describeStore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
collect basic descriptive statistics on ciseStore instances
}
\description{
collect basic descriptive statistics on ciseStore instances
}
\usage{
describeStore(st, genetag = "probeid", snptag = "snp", ids = NULL, resfilter = force, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{st}{
instance of \code{\link{ciseStore-class}}
}
  \item{genetag}{
string for field name for name of quantitatively assayed feature,
defaults to "probeid"; for GTEx application "gene" is used
}
  \item{snptag}{
string for field name for name of genotype feature
}
  \item{ids}{
integerish vector of ids, can be left NULL to survey entire store
}
  \item{resfilter}{
function applied to job results prior to summarization, defaults to force()
}
  \item{\dots}{
not currently used
}
}
\details{
uses parallel infrastructure of BatchJobs on contents managed by st@reg
%%  ~~ If necessary, more details than the description above ~~
}
\value{
list with elements ntests, ngene.uniq, nsnp.uniq
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
%\note{
%%%  ~~further notes~~
%}
%
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
 library(geuvStore)
 mm = makeGeuvStore()
 describeStore(mm, ids=partialIds()[1:10], resfilter=function(x) x[x$mindist < 50000])
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
