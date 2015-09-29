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
Objects can be created by calls of the form \code{new("ciseStore", reg=reg, ...)}.  All arguments must be named.

We can also use
\code{ciseStore(reg, validJobs, addProbeMap = TRUE, addRangeMap = TRUE)}
and the probemap and rangeMap slots will be populated
appropriately.  If validJobs is missing, the validJobs slot will be
populated by findDone(reg).  This may be problematic for handcrafted
extracts from archives.
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

Function \code{describeStore} uses batchMapResults and reduceResults
to leverage a parallel environment to collect information on numbers
of tests and features.  Arguments are described in the associated
man page.
}
%\references{
%%  ~~put references to the literature/web site here~~
%}
%\author{
%%%  ~~who you are~~
%}
\note{
the construction of the maps occurs via \code{\link{storeApply}}, which
%will use \code{\link[BiocParallel]{bplapply}} and thus the currently
%registered bpparam plays a role in determining the speed of construction.
will use \code{\link[foreach]{foreach}}, so that registration
of a parallel back end
using, e.g., \code{\link[doParallel]{registerDoParallel}},
will determine the speed of construction

Any registry job results that do not inherit from GRanges are
mapped to NULL and will not be present in ultimate maps.
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
 store = makeGeuvStore()
 store
 }
}
\keyword{classes}
