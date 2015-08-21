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

storeToFf = function( store, field, ids=NULL, 
   filter=force, ..., checkField=FALSE, ischar = FALSE ) {
#
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
        f <- ffapp2(f, x)
        oldw = options()$warn
        options(warn=-1)
        delete(x)
        rm(x)
        options(warn=oldw)
    }
    f
}

  stopifnot(length(field)==1 && is.character(field))
  if (is.null(ids)) ids = store@validJobs
  if (checkField) {
       result1 = loadAndFilterResult(reg=store@reg, id=ids[1],filter=filter)
       stopifnot(field %in% names(mcols(result1)))
       if (is.character(mcols(result1)[,field]) & !ischar) {
           message("note: checkField identifies entity as character but ischar == FALSE; setting to TRUE")
           ischar = TRUE
           }
       }
##BP  tmp = bplapply(ids, function(x) { #}
      tmp = foreach(x=ids) %dopar% {
        patt = paste0("ff_", x)
        g = as.numeric
        if (ischar) g = function(x) factor(as.character(x))
        ff(g(mcols(loadAndFilterResult(reg=store@reg, id=x, filter=filter))[[field]]), pattern=patt)
##BP      })
        }
  suppressMessages({do.call(cleanc, tmp)})
}

extractByProbes = function(store, probeids, extractTag="probeid") {
  pmap = getProbeMap(store)
  if (any(is.na(probeids))) {
     message("omitting some NA probeids...")
     probeids = as.character(na.omit(probeids))
     }
  uids = unique(pmap[ match(probeids, pmap[,1]), 2 ])
  uids = as.integer( na.omit(uids) )
##BP  ans = bplapply( uids, function(x) {
  ans = foreach (x=uids) %dopar% {
       tmp = getResult(store, x)  # thinner than getResults on all ids
       if (length(tmp)>0) tmp$jobid = x
       tmp[ which(mcols(tmp)[[extractTag]] %in% probeids) ]
##BP       })
       }
  unlist(GRangesList(ans))  # seems a nuisance
}

extractBySymbols = function(store, symbols, sym2probe, ...) {
#
# sym2probe is named vector c(sym1=p1, sym2=p2, and so on)
#
 stopifnot(is(sym2probe, "character"), is(names(sym2probe), "character"))
 rmap = names(sym2probe)
 names(rmap) = as.character(sym2probe)
 ans = extractByProbes(store, sym2probe[symbols], ...)
 if ("sym" %in% names(mcols(ans))) message("clobbering 'sym' element of mcols of result")
 ans$sym =  rmap[ ans$probeid ]
 ans
}
 
extractByRanges = function(store, gr) {
  rmap = getRangeMap(store)
  fi = findOverlaps( rmap, gr )
  sh = queryHits(fi)
  ids = as.integer(unique(rmap[sh]$jobid))
##BP  ans = bplapply(ids, function(x) {
  ans = foreach(x=ids) %dopar% {
      tmp = subsetByOverlaps(getResult(store, x),
       gr)
      if (length(tmp) == 0) return(tmp)
      tmp$jobid = x
      tmp
##BP      })
  }
  ans = ans[ which(sapply(ans,length)>0) ]
  stopifnot(length(ans)>0)
  GenomicRanges::unlist(GRangesList(ans))
}


storeApply = function( store, f, n.chunks, ids=NULL, ... , verbose=FALSE ) {
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
 curids = getJobIds( store )
 if (!is.null(ids)) ids = intersect(ids,curids)
 chs = getStoreIDchunks( store, n.chunks, ids=ids ) #chunk( ids, n.chunks = n.chunks )
# probably need to intersect chs with ids or ids is ignored
 fOnRetrieval = function(ch) reduceResultsList( getRegistry(store), ch,
      fun=function(job, res) f(res) )
##BP bplapply( chs, fOnRetrieval, ... )
 foreach(x=chs) %dopar% fOnRetrieval(x, ...)
}

makeProbeMap = function(store, ..., probetag="probeid") {
 chk1 = loadResult( store@reg, 1)
 stopifnot(probetag %in% names(mcols(chk1)))
 plist = storeApply( store, function(x) unique(as.character(mcols(x)[[probetag]])), ... )
 ul = unlist(plist, recursive=FALSE)
 lens = sapply(ul, length)
 jobn = as.numeric(names(ul))
 jobnum = rep(jobn, lens)
 data.frame(probeid=unlist(ul), jobnum=jobnum, stringsAsFactors=FALSE)
}

makeRangeMap = function(store, ...) {
 chk1 = loadResult( store@reg, 1)
 stopifnot(is(chk1, "GRanges"))
 plist = storeApply( store, range )  # storeApply will create a list of lists
 ul = unlist(GRangesList(unlist(plist)))
 ul$jobid = names(ul)
 ul
}

getStoreIDchunks = function( store, n.chunks, ids=NULL ) {
##BP if (missing(n.chunks)) n.chunks = bpworkers(bpparam())
 if (missing(n.chunks)) n.chunks = getDoParWorkers()
 curids = getJobIds(store)
 if (!is.null(ids)) curids = intersect(ids, curids)
 chunk( curids, n.chunks = n.chunks )
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

DFstoreToFf = function( store, field, ids=NULL, 
   ..., checkField=FALSE, ischar = FALSE ) {
#
# for BySNP assessment, we saved data.frame instances only
#
stopifnot( inherits(store, "Registry") )
cleanc = function (...) 
{
#
# avoids problems with > 1500 ff open
# NOTE: simple modification to ffbase::c.ff 0.11.3, GPL-3
#
    l <- list(...)
    f <- NULL
    for (x in l) {
        f <- ffapp2(f, x)
        delete(x)
        rm(x)
    }
    f
}

  stopifnot(length(field)==1 && is.character(field))
  if (is.null(ids)) ids = store@validJobs
  if (checkField) {
       result1 = loadResult(reg=store, id=ids[1]) #,filter=filter)
       stopifnot(field %in% names(result1))
       if (is.character(result1[,field]) & !ischar) {
           message("note: checkField identifies entity as character but ischar == FALSE; setting to TRUE")
           ischar = TRUE
           }
       }
##BP  tmp = bplapply(ids, function(x) {
     tmp = foreach(x=ids) %dopar% {
      patt = paste0("ff_", x)
      g = as.numeric
      if (ischar) g = function(x) factor(as.character(x))
      #ff(g(loadResult(reg=store, id=x, filter=filter)[[field]]), pattern=patt)
      ff(g(loadResult(reg=store, id=x)[[field]]), pattern=patt)
#BP      })
    }
  suppressMessages({do.call(cleanc, tmp)})
}
