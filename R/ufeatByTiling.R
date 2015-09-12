
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
      tochop = unlist( ans[ bige ] )
      chs = chunk(1:length(tochop), chunk.size=maxlen)
      chopped = lapply(chs, function(x) tochop[x])
      ans = c(okl, chopped)
      }
   ans
}

