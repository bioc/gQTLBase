
describeStore.old = function(st, genetag="probeid", snptag="snp", ...) {
 dd = storeApply(st,function(x)
                            list(jlen=length(x),
                                 ugn = unique(mcols(x)[[genetag]]),
                                 usn = unique(mcols(x)[[snptag]]),
                                 gps = sapply(split(mcols(x)[[genetag]], mcols(x
)[[snptag]]), length),
                                 spg = sapply(split(mcols(x)[[snptag]], mcols(x)
[[genetag]]), length)), flatten1=TRUE, ...)
#
# some of the following are slow.  isn't clear that parallelization
# with typical tools is helpful
#
 ntests = sum(unlist(sapply(dd, function(x) x$jlen)))
 ngene = length(unique(unlist(sapply(dd, function(x) x$ugn))))
 nsnp = length(unique(unlist(sapply(dd, function(x) x$usn))))
 gpshist = hist(unlist(lapply(dd, function(x)x$gps)), plot=FALSE)
 spghist = hist(unlist(lapply(dd, function(x)x$spg)), plot=FALSE)
 list(ntests=ntests, ngene=ngene, nsnp=nsnp, gpshist=gpshist, spghist=spghist)
}

#describeStore2 = function(st  # use BatchMapResults directly?

describeStoreToTmpReg = function(st, genetag="probeid", snptag="snp", ids=NULL,
  resfilter=force, ...) {
  dir.create(tfolder <- tempfile())
  tmpreg = makeRegistry(basename(tempfile()), tfolder, packages="gQTLstats")
# counting and listing in parallel
  descJob = function(genetag, snptag) function(job, res, ...) {
    if (is.na(res)) return(NULL)
    res = resfilter(res)
    list(len = length(res), ugn=unique(mcols(res)[[genetag]]),
         usn=unique(mcols(res)[[snptag]]))
    }
  if (!is.null(ids)) batchMapResults(st@reg, tmpreg, descJob(genetag, snptag), ids=ids)
  else batchMapResults(st@reg, tmpreg, descJob(genetag, snptag))
  submitJobs(tmpreg)
  tmpreg
}
describeStore = function(st, genetag="probeid", snptag="snp", ids=NULL, resfilter=force, ...) {
   tmpreg = describeStoreToTmpReg( st=st, genetag=genetag, snptag=snptag, ids=ids, resfilter=resfilter, ...)
   waitForJobs(tmpreg)
   redr = function(aggr, job, res, ...) {
        aggr$len = aggr$len + res$len
        aggr$ugn = unique(c(aggr$ugn, res$ugn))
        aggr$usn = unique(c(aggr$usn, res$usn))
        aggr
        }
   rr = reduceResults(tmpreg, fun=redr)
   list(ntests=rr$len, ngene.uniq = length(rr$ugn), nsnp.uniq=length(rr$usn))
}
