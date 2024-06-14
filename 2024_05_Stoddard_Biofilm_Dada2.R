### DADA2 PREPARING ASVS AND TAXA TABLES

## Load Libraries
library(dada2)
library(csv)
library(tidyverse)

### First Sequencing Run

## Path to fastq Files
workpath <- "/data/home/htstodda/Biofilm_Fastq/workpath"
filepath <- "/data/home/htstodda/Biofilm_Fastq/Biofilm_Code/Biofilm_Scripts/Final Files/2024_05_Stoddard_Biofilm_Fastq1/Sequencing_Run_1" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(filepath)

## Forward and Reverse fastq Filenames Have Format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(filepath, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(filepath, pattern="_R2_001.fastq.gz", full.names = TRUE))

## Extract Sample Names, Assuming Filenames Have Format: SAMPLENAME_SXXX.fastq.gz
sample.names <- sapply(strsplit(basename(fnFs), "_S"), `[`, 1)
head(sample.names)

anyDuplicated(sample.names)

sample_names <- as.data.frame(sample.names)
head(sample_names)

write_csv(sample_names, "/data/home/htstodda/Biofilm_Fastq/workpath/Biofilm_sample_names.csv") # Save Sample Names to a File

## Visualize Quality Profiles
plotQualityProfile(fnFs[1:2])           #forward reads
plotQualityProfile(fnRs[1:2])           #reverse reads

## Place Filtered Files in "filtered", a Subdirectory
filtFs <- file.path(workpath, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(workpath, "filtered", paste0(sample.names, "_R_filt.fastq"))

## Filtering Parameters: 280, 180
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(280, 180),
                     maxN=0, maxEE=Inf, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, trimLeft = 12) # On Windows set multithread=FALSE
#View(out)

read.stats <- out %>%
  as_tibble() %>%
  summarise(
    max.in = max(reads.in),
    min.in = min(reads.in),
    mean.in = mean(reads.in),
    median.in = median(reads.in),
    sum.in = sum(reads.in),
    count.in = sum(ifelse(reads.in > 0, 1, 0)),
    under1000.in = sum(ifelse(reads.in < 1000, 1, 0)),
    max.out = max(reads.out),
    min.out = min(reads.out),
    mean.out = mean(reads.out),
    median.out = median(reads.out),
    sum.out = sum(reads.out),
    count.out = sum(ifelse(reads.out > 0, 1, 0)),
    under1000.out = sum(ifelse(reads.out < 1000, 1, 0))
  )
#View(t(read.stats))

write.csv(out, "/data/home/htstodda/Biofilm_Fastq/workpath/out.csv") # Save to a File

## Learn About Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

## Plot Errors
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

## Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

## Name the derep-class Objects by the Sample Names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

## Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

## Inspecting dada-class Object
dadaFs[[1]]
dadaRs[[1]]

## Merge Paired Reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

## Inspect the Merger data.frame From the First Sample
head(mergers[[1]])

## Construct Sequence Table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

## Inspect Distribution of Sequence Lengths
table(nchar(getSequences(seqtab)))

## Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

## Assign Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/data/home/sharedDatabases/Silva/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxa # Removing Sequence rownames for Display Only
rownames(taxa.print) <- NULL
head(taxa.print)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

## Save Sequence Table and Taxa Table to Files
write_rds(seqtab.nochim, "/data/home/htstodda/Biofilm_Fastq/workpath/seqtab.rds")     #sequence table
write_rds(taxa.print, "/data/home/htstodda/Biofilm_Fastq/workpath/taxtab.rds")     #taxa table


### Repeat for Second Sequencing Run

## Path to fastq Files
filepathRR <- "/data/home/htstodda/Biofilm_Fastq/Biofilm_Code/Biofilm_Scripts/Final Files/2024_05_Stoddard_Biofilm_Fastq1/Sequencing_Run_2" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(filepathRR)

## Forward and Reverse fastq Filenames Have Format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFsRR <- sort(list.files(filepathRR, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRsRR <- sort(list.files(filepathRR, pattern="_R2_001.fastq.gz", full.names = TRUE))

## Extract Sample Names, Assuming Filenames Have Format: SAMPLENAME_SXXX.fastq.gz
sample.names <- sapply(strsplit(basename(fnFsRR), "_S"), `[`, 1)
head(sample.names)

anyDuplicated(sample.names)

sample_namesRR <- as.data.frame(sample.names)
head(sample_namesRR)

write_csv(sample_namesRR, "/data/home/htstodda/Biofilm_Fastq/workpath/Biofilm_sample_namesRR.csv") # Save Sample Names to a File

## Visualize Quality Profiles
plotQualityProfile(fnFsRR[1:2])           #forward reads
plotQualityProfile(fnRsRR[1:2])           #reverse reads

## Place Filtered Files in "filtered", a Subdirectory
filtFsRR <- file.path(workpath, "filteredRR", paste0(sample.names, "_F_filt.fastq"))
filtRsRR <- file.path(workpath, "filteredRR", paste0(sample.names, "_R_filt.fastq"))

## Filtering Parameters: 280, 180 
outRR <- filterAndTrim(fnFsRR, filtFsRR, fnRsRR, filtRsRR, truncLen = c(280, 180),
                     maxN=0, maxEE=Inf, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, trimLeft = 12) # On Windows set multithread=FALSE
#View(outRR)

read.stats <- outRR %>%
  as_tibble() %>%
  summarise(
    max.in = max(reads.in),
    min.in = min(reads.in),
    mean.in = mean(reads.in),
    median.in = median(reads.in),
    sum.in = sum(reads.in),
    count.in = sum(ifelse(reads.in > 0, 1, 0)),
    under1000.in = sum(ifelse(reads.in < 1000, 1, 0)),
    max.out = max(reads.out),
    min.out = min(reads.out),
    mean.out = mean(reads.out),
    median.out = median(reads.out),
    sum.out = sum(reads.out),
    count.out = sum(ifelse(reads.out > 0, 1, 0)),
    under1000.out = sum(ifelse(reads.out < 1000, 1, 0))
  )
#View(t(read.stats))

write.csv(outRR, "/data/home/htstodda/Biofilm_Fastq/workpath/outRR.csv") # Save to a File

## Learn About Error Rates
#errFRR <- learnErrors(filtFsRR, multithread=TRUE)
#errRRR <- learnErrors(filtRsRR, multithread=TRUE)

## To Improve Error Rates for the Second Sequencing Run, an Alternate Method was Used https://github.com/benjjneb/dada2/issues/1307
loessErrfun_mod4 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # jonalim's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),degree = 1, span = 0.95)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

## Check What the Error Rate Looks Like
errFRR <- learnErrors(
  filtFsRR,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod4,
  verbose = TRUE
)

errRRR <- learnErrors(
  filtRsRR,
  multithread = TRUE,
  nbases = 1e10,
  errorEstimationFunction = loessErrfun_mod4,
  verbose = TRUE
)

## Plot Errors
plotErrors(errFRR, nominalQ=TRUE)
plotErrors(errRRR, nominalQ=TRUE)

## Dereplication
derepFsRR <- derepFastq(filtFsRR, verbose=TRUE)
derepRsRR <- derepFastq(filtRsRR, verbose=TRUE)

## Name the derep-class Objects by the Sample Names
names(derepFsRR) <- sample.names
names(derepRsRR) <- sample.names

## Sample Inference
dadaFsRR <- dada(derepFsRR, err=errFRR, multithread=TRUE)
dadaRsRR <- dada(derepRsRR, err=errRRR, multithread=TRUE)

## Inspecting dada-class Object
dadaFsRR[[1]]
dadaRsRR[[1]]

## Merge Paired Reads
mergersRR <- mergePairs(dadaFsRR, derepFsRR, dadaRsRR, derepRsRR, verbose=TRUE)

## Inspect the Merger data.frame From the First Sample
head(mergersRR[[1]])

## Construct Sequence Table
seqtabRR <- makeSequenceTable(mergersRR)
dim(seqtabRR)

## Inspect Distribution of Sequence Lengths
table(nchar(getSequences(seqtabRR)))

## Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtabRR, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtabRR)

## Assign Taxonomy
taxaRR <- assignTaxonomy(seqtab.nochim, "/data/home/sharedDatabases/Silva/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxaRR # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

getN <- function(x) sum(getUniques(x))
track <- cbind(outRR, sapply(dadaFsRR, getN), sapply(dadaRsRR, getN), sapply(mergersRR, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

## Save Sequence Table and Taxa Table to Files
write_rds(seqtab.nochim, "/data/home/htstodda/Biofilm_Fastq/workpath/seqtabRR.rds")     #sequence table
write_rds(taxa.print, "/data/home/htstodda/Biofilm_Fastq/workpath/taxtabRR.rds")        #taxa table


### Combine the Sequence and Taxa Tables from the Two Runs
## Sequence Tables
s1 <- readRDS("/data/home/htstodda/Biofilm_Fastq/workpath/seqtab.rds")
s2 <- readRDS("/data/home/htstodda/Biofilm_Fastq/workpath/seqtabRR.rds")
s <- mergeSequenceTables(s1, s2)
saveRDS(s, "/data/home/htstodda/Biofilm_Fastq/workpath/seqtabFINAL.rds")

## Taxa Tables
t1 <- readRDS("/data/home/htstodda/Biofilm_Fastq/workpath/taxtab.rds")
t2 <- readRDS("/data/home/htstodda/Biofilm_Fastq/workpath/taxtabRR.rds")
t <- rbind(t1, t2)
saveRDS(t, "/data/home/htstodda/Biofilm_Fastq/workpath/taxtabFINAL.rds")
