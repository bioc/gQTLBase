cb2range = function(cbspec="17q12", annopkg="Homo.sapiens", na.action=na.omit) {
  requireNamespace(annopkg)
  pk = get(annopkg) # package name is object name
  allr = select(pk, keys=cbspec, keytype="MAP", columns=c("TXCHROM", "TXSTART", "TXEND"))
  stopifnot(length(allr)>0)
  allr = na.action(allr)
  tmp = GRanges(allr$TXCHROM, IRanges(abs(allr$TXSTART), abs(allr$TXEND)),
         strand=ifelse(allr$TXSTART>0, "+", "-"))
  ans = range(tmp)
  ans$spec = paste(cbspec, collapse=",")
  ans
}
