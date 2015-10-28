
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
describeStore.slow = function(st, genetag="probeid", snptag="snp", ids=NULL, resfilter=force, ...) {
   tmpreg = describeStoreToTmpReg( st=st, genetag=genetag, snptag=snptag, ids=ids, resfilter=resfilter, ...)
   waitForJobs(tmpreg)
   summs = loadResults(tmpreg) # lengths and modest length vectors of strings
   ntests = sum(sapply(summs, function(x) x$len)) # fast
   n.gene.uniq = length(unique(unlist(lapply(summs, function(x) x$ugn)))) # fast?
   n.snp.uniq = length(unique(unlist(lapply(summs, function(x) x$usn)))) # fast?
#   redr = function(aggr, job, res, ...) {
#        aggr$len = aggr$len + res$len
#        aggr$ugn = unique(c(aggr$ugn, res$ugn))
#        aggr$usn = unique(c(aggr$usn, res$usn))
#        aggr
#        }
#   rr = reduceResults(tmpreg, fun=redr)  # sequential, too slow
#   list(ntests=rr$len, ngene.uniq = length(rr$ugn), nsnp.uniq=length(rr$usn))
    list(ntests=ntests, n.gene.uniq=n.gene.uniq, n.snp.uniq=n.snp.uniq)
}

dendroReduce.bj = function(llike, binfun) {
  n = length(llike)
  indl = BBmisc::chunk(1:n, chunk.size=2)
  tname = tempfile()
  tfold = dir.create(tname)
  tmpr = makeRegistry(basename(tname), file.dir=tfold)
  batchMap(tmpr, binfun, indl)
  submitJobs(tmpr, binfun)
  waitForJobs(tmpr)
  tmpr
}

dendroReduce.fe = function(llike, binfun) {
# binfun must be endomorphic
  n = length(llike)
  if (n==1) return(llike)
  indl = BBmisc::chunk(1:n, chunk.size=2)
  red = foreach(i = 1:length(indl)) %dopar%{ binfun(llike[[indl[[i]][1]]],
           llike[[indl[[i]][2]]])}
  Recall(red, binfun)
}
  
describeByFilts = function( st, filtlist, ... ) {
  ds = lapply(filtlist, function(f) describeStore(st, resfilter=f, ...))
  do.call(rbind, ds)
}

jobsByChrom = function(st) {
 rmap = st@rangeMap
 sn = unique(seqnames(rmap))
 lapply(sn, function(x) as.numeric(rmap[seqnames(rmap)==x]$jobid))
}

n.uniq.snp = function(st, snptag="snp", resfilter=force, ids=NULL) {
 jbc = jobsByChrom(st)
 if (!is.null(ids)) {
    jbc = lapply(jbc, function(x) intersect(ids, x))
    jl = sapply(jbc,length)
    if (all(jl==0)) stop("none of the ids selected are available")
    jbc = jbc[which(jl>0)]
    }
 alls = lapply(1:length(jbc), function(i)
    unique(unlist(storeApply(st, function(x)unique(mcols(resfilter(x))[[snptag]]), ids=jbc[[i]])))
 )
 length(unique(unlist(alls)))
}



.describeStore = function(st, genetag = "probeid", snptag = "snp", ids = NULL, resfilter = force, doChecks=TRUE, ...) {
  chkfun = function(x) {
    bad = c(reqsize=NA, reqsat=0, litenloc=NA, len=0)
    if (!is(x, "GRanges")) return(bad)  # presumably NULL
    m = metadata(x)
    c(reqsize=m$requestSize, reqsat=m$nRequestsSatisfied, litenloc=
         m$dimliteGT[1], len=length(x))
  }
  chkstr = NULL
  if (doChecks) chkstr = storeApply( st, chkfun )
  ntests = sum(unlist(storeApply(st, function(x) length(resfilter(x)), flatten1=TRUE, ids=ids)))
  n.gene.uniq <- length(unique(unlist(storeApply(st, f=function(x) unique(mcols(resfilter(x))[[genetag]]), flatten1=TRUE, ids=ids))))
  n.snp.uniq = n.uniq.snp(st, snptag=snptag, resfilter=resfilter)
  list(basic=c(ntests=ntests, n.gene.uniq=n.gene.uniq, n.snp.uniq=n.snp.uniq),
      checks=chkstr)
}

setClass("storeDescription", representation(
    basic="ANY", reqstat="numeric", lenstat="numeric", full="list",
    reqfail="ANY", locfail="ANY"))
setMethod("show", "storeDescription", function(object){
cat("storeDescription:\n")
print(object@basic)
cat("% requests satisfied: ", 100*round(object@reqstat,3), "\n")
cat("% lengths verified: ", 100*round(object@lenstat,3), "\n")
})

describeStore = function(st, genetag = "probeid", snptag = "snp", ids = NULL, resfilter = force, doChecks=TRUE, ...) {
 d1 = .describeStore(st=st, genetag=genetag, snptag=snptag, ids=ids,
      resfilter=resfilter, doChecks=doChecks, ...)
 if (!is.null(d1$checks)) {
  chks = sapply(unlist(d1$checks, recursive=FALSE), force)
  reqstat = mean(breq <- (chks[1,]==chks[2,]), na.rm=TRUE)
  lenstat = mean(lreq <- (chks[3,]==chks[4,]), na.rm=TRUE)
  return(new("storeDescription", basic=d1$basic, reqstat=reqstat,
   lenstat=lenstat, full=d1, reqfail=which(!breq), locfail=which(!lreq)))
  }
 return(new("storeDescription", basic=d1$basic))
}

