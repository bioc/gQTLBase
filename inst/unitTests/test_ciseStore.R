checkCisEstore = function() {
#
# verify that the binary results data are as expected (one example)
#
 require(digest)
 rd1 = system.file("extractP6pop-files/jobs/01/1-result.RData", package="geuvStore2")
 checkTrue(digest(rd1, file=TRUE) == "0600f76369cc70b58c5fe73acce8fe4b")

#
# setup known properties
#

 require(GenomeInfoDb)
 library(doParallel)
 library(parallel)
 registerDoParallel(cores=(detectCores()-1))

# basic operations for a ciseStore instance

 nn = makeGeuvStore2()
 checkTrue(inherits(nn@reg, "Registry"))
 nnrun = ciseStore(reg=nn@reg, 1:160, TRUE, TRUE)
 checkTrue( seqlevelsStyle(nn@rangeMap) == "NCBI" )
 checkTrue( all(genome(nn@rangeMap) == "hg19") )

# extractByProbes
 myp = c("ENSG00000183814.10", "ENSG00000174827.9")
 ebp = extractByProbes(nn, myp)
 checkTrue( length(ebp) == 8888 )
 checkTrue( sum(is.na(ebp$chisq)) == 0 )

# extractByRanges
 rr = range(ebp)
 ebr = extractByRanges( nn, rr )
 checkTrue( length(ebr) == 230141 )
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
onl = structure(c(47630L, 34192L, 33160L, 35274L, 38420L, 26491L, 25649L, 
40555L, 35397L, 41355L, 41513L, 40856L, 33601L, 30242L, 35918L, 
39493L, 19852L, 34195L, 41098L, 44640L, 36420L, 38143L, 30646L, 
44028L, 40670L, 33485L, 29442L, 51983L, 39648L, 28459L, 40266L, 
43562L, 27929L, 31425L, 45772L, 42175L, 33504L, 34002L, 36496L, 
46119L, 46804L, 42225L, 40702L, 35916L, 23941L, 36114L, 33830L, 
36241L, 42881L, 33045L, 43470L, 39853L, 37298L, 32662L, 37051L, 
40763L, 33674L, 32739L, 46722L, 37678L, 25739L, 25182L, 40577L, 
27081L, 49713L, 47464L, 38608L, 32660L, 49990L, 38660L, 31498L, 
37059L, 50536L, 36459L, 27141L, 27487L, 38593L, 42053L, 28843L, 
37092L, 44258L, 44452L, 40274L, 42967L, 54080L, 47933L, 39107L, 
32180L, 38466L, 48392L, 43552L, 31517L, 30740L, 37995L, 32422L, 
43437L, 44452L, 27784L, 40173L, 35067L, 29520L, 35384L, 42113L, 
39063L, 41934L, 48779L, 47066L, 32812L, 34559L, 38690L, 28131L, 
35320L, 49820L, 25250L, 26022L, 34344L, 31685L, 41366L, 50512L, 
41007L, 51314L, 40342L, 41684L, 42670L, 29015L, 32299L, 43499L, 
37378L, 29826L, 34264L, 34697L, 37313L, 81629L, 37219L, 48555L, 
68692L, 131671L, 30899L, 34952L, 53145L, 32153L, 36946L, 47530L, 
42830L, 27834L, 37935L, 49055L, 35496L, 36360L, 42795L, 33667L, 
33792L, 23932L, 30725L, 36331L, 37435L, 44975L, 28115L, 33911L, 
37932L), .Names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", 
"10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", 
"21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", 
"32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", 
"43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", 
"54", "55", "56", "57", "58", "59", "60", "61", "62", "63", "64", 
"65", "66", "67", "68", "69", "70", "71", "72", "73", "74", "75", 
"76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", 
"87", "88", "89", "90", "91", "92", "93", "94", "95", "96", "97", 
"98", "99", "100", "101", "102", "103", "104", "105", "106", 
"107", "108", "109", "110", "111", "112", "113", "114", "115", 
"116", "117", "118", "119", "120", "121", "122", "123", "124", 
"125", "126", "127", "128", "129", "130", "131", "132", "133", 
"134", "135", "136", "137", "138", "139", "140", "141", "142", 
"143", "144", "145", "146", "147", "148", "149", "150", "151", 
"152", "153", "154", "155", "156", "157", "158", "159", "160"
))

 uu = ll[names(onl)]
 checkTrue(all.equal(uu, onl))

# storeApply with ids
 ll = unlist(storeApply(nn, length, ids=c(1:3,103,160)))
limonl = structure(c(47630L, 34192L, 33160L, 42113L, 37932L), .Names = c("1", 
"2", "3", "103", "160"))
 uu = ll[names(limonl)]
 checkTrue(all.equal(uu, limonl))

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
 checkTrue(length(st) == 114982)

# internal
# set.seed(1234)
# x = list(runif(20), runif(20))
# fff = gQTLBase:::ff_from_list(x)
# checkTrue(length(fff)==40)

}
checkCisEstore()
