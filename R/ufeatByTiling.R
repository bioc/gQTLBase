
ufeatByTiling = function(se, tiling, maxlen=20) {
#
# obtain a list of rownames from SummarizedExperiment
# partitioned by tiling, with no repetitions if rowRanges is cut
# by tile boundaries (feature assigned to tile where it appears first)
#
# added feature -- elements length > maxlen are broken up
# to have that length
#
#   tiling = sort(tiling)
#   se = sort(se)
   fo = findOverlaps(se, tiling)
   rn = rownames(se)
   map = data.frame(qh=rn[queryHits(fo)], sh=subjectHits(fo),
       stringsAsFactors=FALSE)
   drop = which(duplicated(map$qh))
   if (length(drop)>0) map = map[-drop,]
   ans = split(map$qh, map$sh)
   lans = sapply(ans,length)
   bige = which(lans > maxlen)
   if (length(bige)>0) {
      okl = ans[ -bige ]
      tochop =  ans[ bige ] 
      sm = vector("list", length=length(tochop))  # ensure lists are from same seq
      for (i in 1:length(tochop) ) {
          chs = chunk(1:length(tochop[[i]]), chunk.size=maxlen)
          sm[[i]] = lapply(chs, function(x) tochop[[i]][x])
          }
      ans = c(okl, unlist(sm, recursive=FALSE))
      }
   ans
}

balancedFeatList = function(se, maxlen=20) {
 sse = split(se, as.character(seqnames(se)))
 ans = lapply(sse, function(x) {
        chunk(rownames(x), chunk.size=maxlen)
        })
 unlist(ans, recursive=FALSE) # return to simple list
}

