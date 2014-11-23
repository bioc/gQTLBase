\name{ciseStore-class}
\Rdversion{1.1}
\docType{class}
\alias{ciseStore-class}
\alias{ciseStore}

\title{Class \code{"ciseStore"}}
\description{
wrap a BatchJobs registry that manages results of a cis-eQTL search
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ciseStore", reg, ...)}.

We can also use
\code{ciseStore(reg, addProbeMap = TRUE, addRangeMap = TRUE)}
and the probemap and rangeMap slots will be populated
appropriately
}
\section{Slots}{
  \describe{
    \item{\code{reg}:}{Object of class \code{"Registry"} BatchJobs
Registry instance }
    \item{\code{validJobs}:}{Object of class \code{"integer"} vector
of valid job identifiers for the registry }
    \item{\code{probemap}:}{Object of class \code{"data.frame"} a map
from expression probe identifiers to job identifiers where results
for the probe are stored}
    \item{\code{rangeMap}:}{Object of class \code{"GRanges"} a map
from ranges on chromosomes, to job identifiers, in mcols()$jobid}
 }
}

\section{Methods}{
show
}
%\references{
%%  ~~put references to the literature/web site here~~
%}
%\author{
%%%  ~~who you are~~
%}
\note{
the construction of the maps occurs via \code{\link{storeApply}}, which
will use \code{\link[BiocParallel]{bplapply}} and thus the currently
registered bpparam plays a role in determining the speed of construction.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("ciseStore")
# get the global assignment back
require(BatchJobs)
if (require(geuvStore)) {
 reg = partialRegistry()
 store = ciseStore(reg, addProbeMap=TRUE, addRangeMap=FALSE)
 store
 }
}
\keyword{classes}
