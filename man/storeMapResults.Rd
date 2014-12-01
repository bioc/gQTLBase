\name{storeMapResults}
\alias{storeMapResults}
\alias{loadAndFilterResult}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
use batchMapResults infrastructure to process
results in a ciseStore instance
}
\description{
use batchMapResults infrastructure to process
results in a ciseStore instance
}
\usage{
storeMapResults(store, reg2, fun, ..., 
    ids = NULL, part = NA_character_, more.args = list())
loadAndFilterResult(reg, 
 id, filter=force, part = NA_character_, missing.ok = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{store}{
an instance of \code{\link[gQTLBase]{ciseStore-class}}
}
  \item{reg}{instance of BatchJobs Registry class}
  \item{reg2}{
an empty instance of the Registry class
(see \code{\link[BatchJobs]{makeRegistry}})
}
  \item{fun}{
A function to map over results in \code{store}, with
formals \code{(job, res, ...)}.
}
  \item{filter}{a function that accepts and returns a GRanges
instance, to be applied just after loading a result from the store}
  \item{\dots}{
additional arguments to vectorize over (should be same length as
\code{length(findDone(store@reg))}
}
  \item{ids}{
ids of job results to be mapped; if missing, map all job results
}
  \item{id}{a single job id}
  \item{part}{
see \code{\link[BatchJobs]{batchMapResults}}
}
  \item{missing.ok}{see \code{\link[BatchJobs]{loadResult}}}
  \item{more.args}{
a list of other arguments to be passed to \code{fun}; default
is empty list.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
integer vector with job ids.  Main purpose is to 
prepare the registry for submitJobs.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
\note{
\code{loadAndFilterResult} is not intended to be exported
and may be removed in future versions.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
if (require(geuvStore)) {
 require(BatchJobs)
 reg = partialRegistry()
 store = ciseStore(reg, addProbeMap=FALSE, addRangeMap=FALSE)
 fd = tempfile()
 tempreg = makeRegistry("tempSMR", file.dir=fd)
 storeMapResults( store, tempreg, fun=function(job, res, ...) length(res) )
 showStatus(tempreg)
 submitJobs(tempreg, 1:2)
 loadResults(tempreg)
 unlink(fd)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
