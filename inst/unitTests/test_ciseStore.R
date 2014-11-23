checkCisEstore = function() {
#
# verify that the binary results data are as expected (one example)
#
  require(digest)
  rd1 = system.file("geuvStore/jobs/01/1-result.RData", package="geuvStore")
  checkTrue(digest(rd1, file=TRUE) == "71a165c73975bbc508ca595002e6c11b")

#
# setup known properties
#

 require(GenomeInfoDb)
 vj = c(1L, 1001L, 101L, 1101L, 1201L, 1301L, 1401L, 1501L, 1601L, 
   1701L, 1801L, 1901L, 2001L, 201L, 2101L, 2201L, 301L, 401L, 501L, 
   601L, 701L, 801L, 901L, 1002L, 102L, 1102L, 1202L, 1302L, 1402L, 
   1502L, 1602L, 1702L, 1802L, 1902L, 2L, 2002L, 202L, 2102L, 2202L, 
   302L, 402L, 502L, 602L, 702L, 802L, 902L, 1003L, 103L, 1103L, 
   1203L, 1303L, 1403L, 1503L, 1603L, 1703L, 1803L, 1903L, 2003L, 
   203L, 2103L, 2203L, 3L, 303L, 403L, 503L, 603L, 703L, 803L, 903L
   )

onl = structure(c(64360L, 74828L, 81351L, 85256L, 80605L, 70500L, 69074L, 
78765L, 76959L, 87299L, 76926L, 75267L, 77434L, 80024L, 74986L, 
85656L, 81924L, 86033L, 91862L, 77926L, 84417L, 78808L, 78878L, 
86451L, 84012L, 87925L, 61949L, 81075L, 72061L, 111395L, 82449L, 
85704L, 67944L, 78054L, 73207L, 76351L, 72877L, 86649L, 71307L, 
75323L, 66124L, 95298L, 94914L, 94346L, 77526L, 93218L, 95213L, 
78858L, 74746L, 75622L, 74128L, 80194L, 76773L, 75153L, 70289L, 
77024L, 73854L, 97933L, 94175L, 81240L, 68090L, 86248L, 71618L, 
104867L, 86626L, 168069L, 78387L, 91499L, 80387L), .Names = c("1", 
"1001", "1002", "1003", "101", "102", "103", "1101", "1102", 
"1103", "1201", "1202", "1203", "1301", "1302", "1303", "1401", 
"1402", "1403", "1501", "1502", "1503", "1601", "1602", "1603", 
"1701", "1702", "1703", "1801", "1802", "1803", "1901", "1902", 
"1903", "2", "2001", "2002", "2003", "201", "202", "203", "2101", 
"2102", "2103", "2201", "2202", "2203", "3", "301", "302", "303", 
"401", "402", "403", "501", "502", "503", "601", "602", "603", 
"701", "702", "703", "801", "802", "803", "901", "902", "903"
))

# basic operations for a ciseStore instance

 reg = geuvStore::partialRegistry()
 checkTrue(inherits(reg, "Registry"))
 nn = ciseStore(reg, TRUE, TRUE)
 checkTrue( all.equal(dim(nn@probemap), c(690,2)) )
 checkTrue( length(nn@rangeMap) == 69 )
 checkTrue( seqlevelsStyle(nn@rangeMap) == "NCBI" )
 checkTrue( all(genome(nn@rangeMap) == "hg19") )
 checkTrue( length(setdiff(nn@validJobs , vj)) == 0 )

# extractByProbes
 myp = c("ENSG00000183814.10", "ENSG00000174827.9")
 ebp = extractByProbes(nn, myp)
 checkTrue( length(ebp) == 10896 )
 checkTrue( sum(is.na(ebp$chisq)) == 0 )

# extractByRanges
 rr = range(ebp)
 ebr = extractByRanges( nn, rr )
 checkTrue( length(ebr) == 190912 )
 checkTrue(length(unique(ebr$jobid)) == 8)

# SnpMatrix4GRanges
#  seems irrelevant for this package?
# require(VariantAnnotation)
# vcfl = system.file("extdata/chr7-sub.vcf.gz", package="VariantAnnotation")
# r1 = readVcf(vcfl, genome="hg19")
# rd = rowData(r1)
# sm = SnpMatrix4GRanges( rd, TabixFile(vcfl) )
# checkTrue(all.equal(dim(sm), c(2,988)))

# storeApply
 ll = unlist(storeApply(nn, length))
 uu = ll[order(names(ll))]
 checkTrue(all.equal(uu, onl))

# storeMapResults
 fd = tempfile()
 tempreg=makeRegistry("tempSMR",file.dir=fd)
 storeMapResults(nn, tempreg, fun=function(job,res,...)length(res))
 submitJobs(tempreg, 1:2)
 waitForJobs(tempreg, timeout=30L)
 checkTrue(length(findNotDone(tempreg))==67)
# checkTrue(all.equal(as.numeric(unlist(loadResults(tempreg))), c(74828, 80605)))
 checkTrue(all.equal(as.numeric(unlist(loadResults(tempreg))), c(64360, 74828)))

# storeToFf
 st = storeToFf( nn, "chisq", ids=1:3 )
 checkTrue(length(st) == 216425)

# internal
 set.seed(1234)
 x = list(runif(20), runif(20))
 fff = gQTLBase:::ff_from_list(x)
 checkTrue(length(fff)==40)

}
