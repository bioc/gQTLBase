
setOldClass("Registry")

setClass("ciseStore", representation(reg="Registry", 
  validJobs="integer", probemap="data.frame", rangeMap="GRanges"))

#setMethod("initialize", "ciseStore", function(.Object, reg, ...) {
#  .Object@reg = reg
#  .Object@probemap = data.frame()
#  .Object@rangeMap = GRanges()
#  .Object@validJobs = getResultIDs(.Object)
#  .Object
#})

#setMethod("initialize", c("ciseStore",
#      "numeric", "data.frame", "GRanges"), function(.Object, reg, vj, df, gr) {
#  .Object@reg = reg
#  .Object@probemap = df
#  .Object@rangeMap = gr
#  .Object@validJobs = vj
#  .Object
#})

ciseStore = function(reg, validJobs, addProbeMap=TRUE, addRangeMap=TRUE, probetag="probeid") {
 if (missing(validJobs)) validJobs = findDone(reg)
 tmp = new("ciseStore", validJobs=validJobs,
       reg=reg, probemap=data.frame(), rangeMap=GRanges())
 if (addProbeMap) {
  message("building probe:job map...")
  tmp@probemap = makeProbeMap(tmp, probetag=probetag)
  message("done.")
  }
 if (addRangeMap) {
  message("building range:job map...")
  tmp@rangeMap = makeRangeMap(tmp)
  message("done.")
 }
 tmp
}

setGeneric("getProbeMap", function(x) standardGeneric("getProbeMap"))
setMethod("getProbeMap", "ciseStore", function(x) {
  if (nrow(x@probemap) == 0) stop("empty probe map, please add one")
  x@probemap
  })
setGeneric("getRangeMap", function(x) standardGeneric("getRangeMap"))
setMethod("getRangeMap", "ciseStore", function(x) {
  if (length(x@rangeMap) == 0) stop("empty range map, please add one")
  x@rangeMap
  })

setGeneric("getRegistry", function(x) standardGeneric("getRegistry"))
setMethod("getRegistry", "ciseStore", function(x) x@reg)

setMethod("show", "ciseStore", function(object) {
 cat("ciseStore instance with", length(object@validJobs),
 "completed jobs.\n")
 if (length(object@validJobs)>0) {
  cat("excerpt from job ", object@validJobs[1],":\n")
  suppressMessages({
  show(loadResult(getRegistry(object), (object@validJobs)[1])[1])
  })
 }
# cat("Table of classes of job results:\n")
# table(unlist(bplapply(object@validJobs, function(x) class(loadResults(object@reg, x)[[1]]))))
})

getJobIds = function(store) store@validJobs

