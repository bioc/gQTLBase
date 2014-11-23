
setOldClass("Registry")

setClass("ciseStore", representation(reg="Registry", 
  validJobs="integer", probemap="data.frame", rangeMap="GRanges"))
setMethod("initialize", "ciseStore", function(.Object, reg, ...) {
  .Object@reg = reg
  .Object@probemap = data.frame()
  .Object@rangeMap = GRanges()
  .Object@validJobs = getResultIDs(.Object)
  .Object
})

ciseStore = function(reg, addProbeMap=TRUE, addRangeMap=TRUE) {
 tmp = new("ciseStore", reg)
 if (addProbeMap) {
  message("building probe:job map...")
  tmp@probemap = makeProbeMap(tmp)
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
setMethod("getProbeMap", "ciseStore", function(x)
  x@probemap)
setGeneric("getRangeMap", function(x) standardGeneric("getRangeMap"))
setMethod("getRangeMap", "ciseStore", function(x)
  x@rangeMap)

setGeneric("getRegistry", function(x) standardGeneric("getRegistry"))
setMethod("getRegistry", "ciseStore", function(x) x@reg)

setMethod("show", "ciseStore", function(object) {
 cat("ciseStore instance with", length(object@validJobs),
 "completed jobs.\n")
# cat("first job result:\n")
# show(loadResults(getRegistry(object), (object@validJobs)[1]))
# cat("Table of classes of job results:\n")
# table(unlist(bplapply(object@validJobs, function(x) class(loadResults(object@reg, x)[[1]]))))
})

getJobIds = function(store) store@validJobs

