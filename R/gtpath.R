
gtpath1kg = function(chrtok) {
gsub("%%N%%", chrtok, 
"http://1000genomes.s3.amazonaws.com/release/20110521/ALL.chr%%N%%.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz")
}
