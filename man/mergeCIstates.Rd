\name{mergeCIstates}
\alias{mergeCIstates}
\alias{mergeGWhits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
merge ChromImpute chromatin states, or GWAS hit
indicators, to a GRanges
}
\description{
merge ChromImpute chromatin states, or GWAS hit
indicators, to a GRanges
}
\usage{
mergeCIstates(gr, ermaset, epig, genome = "hg19", useErma = TRUE)
mergeGWhits(gr, gwcat, use=c("both", "addr", "name")[1],
     grSnpField="SNP")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gr}{
a GRanges instance
}
  \item{ermaset}{
an instance of \code{\link[erma]{ErmaSet-class}}
}
  \item{gwcat}{
an instance of \code{\link[gwascat]{gwaswloc-class}}, or any
compliant GRanges instance -- must have mcols field \code{SNPS}
with snp identifier
}
  \item{epig}{
the standardized epigenome name of the epigenome to use
}
  \item{genome}{
a tag for genome build
}
  \item{useErma}{
logical -- at the moment, must be TRUE; plan is to allow use
of elements of AnnotationHub
}
  \item{use}{
character string selecting approach for linking loci in \code{gr}
to those in \code{gwcat} -- if \code{"both"}, coincidence in address
or name are both checked and used; if \code{"addr"}, only address
is checked, if \code{"name"}, only SNP name.
}
  \item{grSnpField}{
character string naming the field in mcols(gr) with SNP id
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
for mergeCIstates,
a GRanges instance with additional fields in mcols: fullStates, states,
and statecols, denoting respectively the full annotation of ChromImpute
for the inferred state, an abbreviated tag that collapses related states,
and a color tag for rendering, that does not replicate the colors in
the ChromImpute bed files.  The states field is a factor with levels
\code{c("Het", "DNAse", "Enh", "Prom", "Quies", "ReprPC", "Tss", "Tx", 
"ZNF/Rp")}.

for  mergeGWhits, a single mcols field is added, \code{isGwasHit}, that
is 1 for coincident hit and 0 otherwise.  Eventually phenotype information
will be collected and added.
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
if (require(gwascat) && require(erma)) {
#
# demonstrate Tx state for exon starts
#
  gm = resize(genemodel("ORMDL3"),1)
  es = makeErmaSet()
  g1 = mergeCIstates(gm, es, "LNG.FET")
  g1
#
# set up for GWAS
#
  require(GenomeInfoDb)
  data(ebicat37)
  genome(ebicat37) = "hg19"
  seqlevelsStyle(ebicat37) = "UCSC"
  g1 = c(g1, g1[1])  # add a known hit
  start(g1[length(g1)]) = 38062196
  mergeGWhits(g1, ebicat37)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
