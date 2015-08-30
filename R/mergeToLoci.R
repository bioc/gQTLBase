
mergeCIstates = function(gr, ermaset = NULL, epig, genome="hg19", importFull=FALSE, useErma=TRUE, stateGR=NULL) {
#
# label each range in a GenomicRanges with the chromatin state
# for a selected epigenome
#
      if (!useErma) stop("only set up for erma; planning AnnotationHub for 2016")
      if (is.null(stateGR)) {
        cd = colData(ermaset)
        ind = match(epig, cd$Epigenome.Mnemonic)
        fn = files(ermaset)[ind]
        if (importFull) st = import(fn, genome=genome)
        else st = import(fn, which=gr, genome=genome)
        }
      else st = stateGR
      ov = findOverlaps(gr, st)
      gr$fullStates[queryHits(ov)] = st$name[subjectHits(ov)]
      abbCIstates = get(load(system.file("data/abbCIstates.rda", package="erma")))
      abbCIcols = get(load(system.file("data/abbCIcols.rda", package="erma")))
stlev = c("Het", "DNAse", "Enh", "Prom", "Quies", "ReprPC", "Tss", "Tx", 
"ZNF/Rp")
      gr$states = factor(abbCIstates[gr$fullStates], levels=stlev)
      gr$statecols = abbCIcols[gr$states]
      gr
  }

mergeGWhits = function(gr, gwcat, use=c("both", "addr", "name")[1],
     grSnpField="SNP") {
 stopifnot(genome(gr)[1] == genome(gwcat)[1])
 # idea is to put a string with gwas hit phenotype in gwcat on coincident loci in gr
 if ("isGwasHit" %in% names(mcols(gr))) warning("will overwrite mcols field 'isGwasHit'")
 gr$isGwasHit = 0
 if (use == "both" | use == "addr") {
   fo = findOverlaps(gr, gwcat)
   gr$isGwasHit[queryHits(fo)] = 1
   }
 if (use == "both" | use == "name") {
   m = na.omit(match(mcols(gr)[[grSnpField]], gwcat$SNPS))
   gr$isGwasHit[m] = 1
   }
 gr
}
