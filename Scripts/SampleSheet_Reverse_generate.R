##!/usr/local/bin/R
##!Rscript to generate the nextflow r 2.0 input sample sheet with reverse strand
## Xiaohui Zhao (xz289@cam.ac.uk, 20/04/2021)

fastq_1.1 <- list.files(path="/rds-d1/project/rsh46/rds-rsh46-ctr-bfx/CTR-Projects/CTR_myt25_0006/NGS_M_Sheridan_40616_run1/",
                      pattern = "*.fastq.gz", recursive=T)
fastq_1.2 <- list.files(path="/rds-d1/project/rsh46/rds-rsh46-ctr-bfx/CTR-Projects/CTR_myt25_0006/NGS_M_Sheridan_40616_run2/",
                        pattern = "*.fastq.gz", recursive=T)
fastq_1.1 <- fastq_1.1[-grep("work|results", fastq_1.1)]
fastq_1.2 <- fastq_1.2[-grep("work|results", fastq_1.2)]

fastq_1.ori   <- c(fastq_1.1,fastq_1.2)

ngroup    <- gsub(".._ds.*/", "", fastq_1.ori)
ngroup    <- gsub("_L0..*", "", ngroup)

group <- paste0(ngroup, "_L", rep(c(1:4), length=272))
replicate <- rep(c(1:2), each =136)
fastq_2 <- rep(" ", length=272)
strandedness <- "reverse"

fastq_1 <- c(paste0("NGS_M_Sheridan_40616_run1/", fastq_1.1),
  paste0("NGS_M_Sheridan_40616_run2/", fastq_1.2))
sampleSheet <- data.frame(group=group, replicate=replicate, fastq_1=fastq_1,
                          fastq_2=fastq_2,strandedness=strandedness)
write.csv(sampleSheet, file = "/rds-d1/project/rsh46/rds-rsh46-ctr-bfx/CTR-Projects/CTR_myt25_0006/SampleSheet_nextflow_reverse.csv", quote=F, row.names=F)



##----------FIN--------------------------------##