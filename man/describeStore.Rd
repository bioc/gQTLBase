\name{describeStore}
\alias{describeStore}
\alias{describeByFilts}
\alias{storeDescription-class}
\alias{show}
\alias{show,storeDescription-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
collect basic descriptive statistics on ciseStore instances
}
\description{
collect basic descriptive statistics on ciseStore instances
}
\usage{
describeStore(st, genetag = "probeid", snptag = "snp", ids = NULL,
                 resfilter = force, doChecks = TRUE, ...)
describeByFilts(st, filtlist, ...)
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
  \item{filtlist}{
a list of functions suitable as resfilter arguments
}
  \item{doChecks}{logical -- if true, will collect information on
match between number probes requested and number reported on,
and two scans of VCF loci in cis to probes.  See details.}
  \item{\dots}{
used with describeByFilts, pass to storeApply
}
}
\details{
uses parallel infrastructure of foreach on contents managed by st@reg

describeByFilts returns a matrix of descriptions with one row
per filtlist element

storeDescription holds results of a describe task and includes
information on noncongruence of features with cis tests
and of results of two distinct scans of VCF: one with readGT
on a single sample, the other with readVcf on all samples.
If there are discrepancies between features given and tests
returned, [storeDescription]@reqfail will give the job ids for these.
If there are discrepancies between the numbers of loci retrieved
on the two VCF scans, @locfail will tive the job ids for these.
@reqfail events may be legitimate when a feature has no SNP in
cis at the given radius.  @locfail events usually indicate an
I/O problem and the jobs should be resubmitted.

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
 library(geuvStore2)
 mm = makeGeuvStore2()
 describeStore(mm, ids=1:10, resfilter=function(x) x[x$mindist < 50000])
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
