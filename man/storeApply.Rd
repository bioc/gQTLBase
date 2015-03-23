\name{storeApply}
\alias{storeApply}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
apply a function
over job results in a ciseStore instance
}
\description{
apply a function 
over job results in a ciseStore instance
}
\usage{
storeApply(store, f, n.chunks, ids=NULL, ..., verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{store}{
instance of \code{\link{ciseStore-class}}
}
  \item{f}{
function on GRanges stored in ciseStore
}
  \item{n.chunks}{
Number of chunks into which the jobs are
to be broken; the
%series of chunks is handed to \code{\link[BiocParallel]{bplapply}}
%and \code{\link[BatchJobs]{reduceResults}} is used to extract
series of chunks is handed to \code{\link[foreach]{foreach}}
to extract
results and apply \code{f} to them.

%If missing, the value of \code{bpworkers(bpparam())} will be
If missing, the value of \code{getDoParWorkers()}
used.
}
  \item{ids}{defaults to NULL; if non-null, the jobs to
be processed are limited to those identified in this vector.}
  \item{\dots}{
additional arguments to foreach
}
  \item{verbose}{
if TRUE will allow progressbars and other messages to display 
}
}
\details{
The chunking of job identifiers will determine the degree
of parallelization of application, and the form of the list that
is returned.
}
\value{
A list whose structure depends on the chunking of job identifiers.
See the examples.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%%% ~put references to the literature/web site here ~
%}
%\author{
%%%  ~~who you are~~
%}
\note{
eqtlStore imports BiocParallel's bpparam function, and this
determines in real time the
number of workers to be employed by storeApply.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{storeMapResults}} will apply over the
store using the batch jobs
submission infrastructure and can
target specific results via \code{ids}; \code{storeApply} uses bplapply
over the entire store
}
\examples{
if (require(geuvStore)) {
 require(BatchJobs)
 reg = partialRegistry()
 store = ciseStore(reg, addProbeMap=FALSE, addRangeMap=FALSE)
 storeApply(store, length)
 storeApply(store, length, ids=c(1:3,603))
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
