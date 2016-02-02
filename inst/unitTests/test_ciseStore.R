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
 library(doParallel)
 library(parallel)
 registerDoParallel(cores=(detectCores()-1))

vj = c(1L, 1001L, 101L, 1101L, 1201L, 1301L, 1401L, 1501L, 1601L,
   1701L, 1801L, 1901L, 2001L, 201L, 2101L, 2201L, 301L, 401L, 501L,
   601L, 701L, 801L, 901L, 1002L, 102L, 1102L, 1202L, 1302L, 1402L,
   1502L, 1602L, 1702L, 1802L, 1902L, 2L, 2002L, 202L, 2102L, 2202L,
   302L, 402L, 502L, 602L, 702L, 802L, 902L, 1003L, 103L, 1103L,
   1203L, 1303L, 1403L, 1503L, 1603L, 1703L, 1803L, 1903L, 2003L,
   203L, 2103L, 2203L, 3L, 303L, 403L, 503L, 603L, 703L, 803L, 903L,
1029L,  129L,   1529L,  1829L,  2129L,  29L,    529L,   829L,
1129L,  1329L,  1629L,  1929L,  2229L,  329L,   629L,   929L,
1229L,  1429L,  1729L,  2029L,  229L,   429L,   729L)


onl = structure(c(64360L, 74828L, 80605L, 78765L, 76926L, 80024L, 81924L, 
77926L, 78878L, 87925L, 72061L, 85704L, 76351L, 71307L, 95298L, 
77526L, 74746L, 80194L, 70289L, 97933L, 68090L, 104867L, 78387L, 
81351L, 70500L, 76959L, 75267L, 74986L, 86033L, 84417L, 86451L, 
61949L, 111395L, 67944L, 73207L, 72877L, 75323L, 94914L, 93218L, 
75622L, 76773L, 77024L, 94175L, 86248L, 86626L, 91499L, 85256L, 
69074L, 87299L, 77434L, 85656L, 91862L, 78808L, 84012L, 81075L, 
82449L, 78054L, 86649L, 66124L, 94346L, 95213L, 78858L, 74128L, 
75153L, 73854L, 81240L, 71618L, 168069L, 80387L, 72875L, 67546L, 
80230L, 70001L, 84541L, 69049L, 86293L, 82150L, 69967L, 84065L, 
74547L, 92132L, 89845L, 83086L, 78478L, 81637L, 73130L, 76907L, 
87246L, 87733L, 84141L, 139669L, 88165L), .Names = c("1", "1001", 
"101", "1101", "1201", "1301", "1401", "1501", "1601", "1701", 
"1801", "1901", "2001", "201", "2101", "2201", "301", "401", 
"501", "601", "701", "801", "901", "1002", "102", "1102", "1202", 
"1302", "1402", "1502", "1602", "1702", "1802", "1902", "2", 
"2002", "202", "2102", "2202", "302", "402", "502", "602", "702", 
"802", "902", "1003", "103", "1103", "1203", "1303", "1403", 
"1503", "1603", "1703", "1803", "1903", "2003", "203", "2103", 
"2203", "3", "303", "403", "503", "603", "703", "803", "903", 
"1029", "1129", "1229", "129", "1329", "1429", "1529", "1629", 
"1729", "1829", "1929", "2029", "2129", "2229", "229", "29", 
"329", "429", "529", "629", "729", "829", "929"))

# basic operations for a ciseStore instance

 reg = geuvStore::partialRegistry()
 checkTrue(inherits(reg, "Registry"))
 nn = ciseStore(reg=reg, validJobs=geuvStore::partialIds(), TRUE, TRUE)
# checkTrue( all.equal(dim(nn@probemap), c(920,2)) )
# checkTrue( length(nn@rangeMap) == 92 )
 checkTrue( "NCBI" %in% seqlevelsStyle(nn@rangeMap) )
 checkTrue( all(genome(nn@rangeMap) == "hg19") )
# checkTrue( length(setdiff(nn@validJobs , vj)) == 0 )

# extractByProbes
 myp = c("ENSG00000183814.10", "ENSG00000174827.9")
 ebp = extractByProbes(nn, myp)
# checkTrue( length(ebp) == 10896 )
 checkTrue( sum(is.na(ebp$chisq)) == 0 )

# extractByRanges
 rr = range(ebp)
 ebr = extractByRanges( nn, rr )
# checkTrue( length(ebr) == 278166 )
# checkTrue(length(unique(ebr$jobid)) == 11)

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
 uu = ll[names(onl)]
# checkTrue(all.equal(uu, onl))

# storeApply with ids
 ll = unlist(storeApply(nn, length, ids=c(1:3,603,903)))
 limonl = structure(c(64360L, 73207L, 78858L, 81240L, 80387L), .Names = c("1", 
"2", "3", "603", "903"))
 uu = ll[names(limonl)]
# checkTrue(all.equal(uu, limonl))

# storeMapResults  -- inducing timeouts, august 2015
# fd = tempfile()
# library(BatchJobs)
# tempreg=makeRegistry("tempSMR",file.dir=fd)
# storeMapResults(nn, tempreg, fun=function(job,res,...)length(res))
# submitJobs(tempreg, 1:2)
# waitForJobs(tempreg, timeout=30L)
# checkTrue(length(findNotDone(tempreg))==90)
# checkTrue(all.equal(as.numeric(unlist(loadResults(tempreg))), c(74828, 80605)))
# checkTrue(all.equal(as.numeric(unlist(loadResults(tempreg))), c(64360, 74828)))

# storeToFf
 st = storeToFf( nn, "chisq", ids=1:3 )
# checkTrue(length(st) == 216425)

# internal
# set.seed(1234)
# x = list(runif(20), runif(20))
# fff = gQTLBase:::ff_from_list(x)
# checkTrue(length(fff)==40)

}
checkCisEstore()
