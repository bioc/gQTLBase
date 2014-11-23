getResultIDs = function(store) {
#
# purpose is to obtain ids that are valid for a registry that
# has had some jobs removed ... such as geuvStore
#
  reg = store@reg
  fdir = reg$file.dir
  jdirs = dir(paste0(fdir, "/jobs"), full.names=TRUE)
  jn = lapply(jdirs, function(x) dir(x, pattern="RData$"))
  as.integer(as.numeric(unlist(lapply(jn, function(x) sub("-.*", "", x)))))
  }

loadAndFilterResult = function(reg, id, filter=force, part = NA_character_, missing.ok = FALSE) {
 filter(loadResult(reg=reg,id=id,part=part,missing.ok=missing.ok))
 }

storeMapResults = function( store, reg2, fun, ..., ids = NULL,
   part = NA_character_, more.args = list() ) {
#
# purpose is to use batchJobs infrastructure to apply a function
# over all results of an existing batch submission
#
# formals of fun must include job, res
#
 stopifnot(inherits(reg2, "Registry"))
 validids = store@validJobs
 reg1 = getRegistry(store)
 if (!is.null(ids)) validids = intersect(validids, ids)
 batchMapResults( reg1, reg2, fun, ids=validids, ..., part=part, more.args=more.args )
}

makeTempRegistry = function(...) {
 tf = tempfile()
 ans = makeRegistry( basename(tf), file.dir=tf, ... )
 Sys.sleep(1L)
 ans
}

killTempRegistry = function(reg) {
 unlink(reg$file.dir, recursive=TRUE)
}

storeToFf = function( store, field, ids=NULL, filter=force, ..., checkField=FALSE ) {
#
# getter must return a numeric vector
#
cleanc = function (...) 
{
#
# avoids problems with > 1500 ff open
# NOTE: simple modification to ffbase::c.ff 0.11.3, GPL-3
#
    l <- list(...)
    f <- NULL
    for (x in l) {
        f <- ffappend(f, x)
        delete(x)
        rm(x)
    }
    f
}

  stopifnot(length(field)==1 && is.character(field))
  if (is.null(ids)) ids = store@validJobs
  if (checkField) {
       result1 = loadAndFilterResult(reg=store@reg, id=ids[1],filter=filter)
       stopifnot(field %in% names(mcols(result1)))
       }
  tmp = bplapply(ids, function(x) {
      patt = paste0("ff_", x)
      ff(as.numeric(mcols(loadAndFilterResult(reg=store@reg, id=x, filter=filter))[[field]]), pattern=patt)
      })
  suppressMessages({do.call(cleanc, tmp)})
}

extractByProbes = function(store, probeids) {
  pmap = getProbeMap(store)
  uids = unique(pmap[ match(probeids, pmap[,1]), 2 ])
  ans = bplapply( uids, function(x) {
       tmp = getResult(store, x)  # thinner than getResults on all ids
       if (length(tmp)>0) tmp$jobid = x
       tmp[ which(tmp$probeid %in% probeids) ]
       })
  unlist(GRangesList(ans))  # seems a nuisance
}
  
 
extractByRanges = function(store, gr) {
  rmap = getRangeMap(store)
  fi = findOverlaps( rmap, gr )
  sh = queryHits(fi)
  ids = as.integer(unique(rmap[sh]$jobid))
  ans = bplapply(ids, function(x) {
      tmp = subsetByOverlaps(getResult(store, x),
       gr)
      if (length(tmp) == 0) return(tmp)
      tmp$jobid = x
      tmp
      })
  ans = ans[ which(sapply(ans,length)>0) ]
  stopifnot(length(ans)>0)
  GenomicRanges::unlist(GRangesList(ans))
}


storeApply = function( store, f, n.chunks, ... , verbose=FALSE ) {
 oldPB = getOption("BBmisc.ProgressBar.style")
 oldBJV = getOption("BatchJobs.verbose")
 on.exit( {
   options("BBmisc.ProgressBar.style"=oldPB)
   options("BatchJobs.verbose"=oldBJV)
   } )
 if (!verbose)  {
   options(BBmisc.ProgressBar.style="off")
   options(BatchJobs.verbose=FALSE)
   }
 ids = getJobIds( store )
 chs = getStoreIDchunks( store, n.chunks ) #chunk( ids, n.chunks = n.chunks )
# probably need to intersect chs with ids or ids is ignored
 fOnRetrieval = function(ch) reduceResultsList( getRegistry(store), ch,
      fun=function(job, res) f(res) )
 bplapply( chs, fOnRetrieval, ... )
}

makeProbeMap = function(store, ...) {
 plist = storeApply( store, function(x) unique(as.character(mcols(x)$probeid)), ... )
 ul = unlist(plist, recursive=FALSE)
 lens = sapply(ul, length)
 jobn = as.numeric(names(ul))
 jobnum = rep(jobn, lens)
 data.frame(probeid=unlist(ul), jobnum=jobnum, stringsAsFactors=FALSE)
}

makeRangeMap = function(store, ...) {
 plist = storeApply( store, range )
 ul = unlist(GRangesList(unlist(plist)))
 ul$jobid = names(ul)
 ul
}

getStoreIDchunks = function( store, n.chunks ) {
 if (missing(n.chunks)) n.chunks = bpworkers(bpparam())
 chunk( getJobIds( store ), n.chunks = n.chunks )
}

getResult = function(store, ind) {
 stopifnot(length(ind)==1)
 stopifnot(ind %in% getJobIds(store) )
 loadResult( store@reg, ind )
}

getResultList = function(store, inds) {
 stopifnot(all(inds %in% getJobIds(store) ))
 loadResults( store@reg, inds )
}

