library(here)
library(Biostrings)
library(ShortRead)
library(dada2); packageVersion("dada2")

setwd("~/project_data/")
path <- "221024_UNC2X_KKYJ5-KKYGJ-NETTI-in/221024_UNC2X_KKYJ5-KKYGJ-NETTI"

list.files(path)
fnF <- list.files(path, pattern="_R1_", full.names = TRUE)
fnR <- list.files(path, pattern="_R2_", full.names = TRUE)

sam <- sapply(strsplit(basename(fnF), "_"), `[`, 1)
if(!identical(sam, sapply(strsplit(basename(fnR), "_"), `[`, 1))) stop("F/R fastqs mismatch") 

#plotQualityProfile(fnF[c(1,10,100)])
#plotQualityProfile(fnR[c(1,10,100)])

plotComplexity(fnF[[13]])
plotComplexity(fnR[[53]])
# Minor "lower" complexity mode, around 13 vs. 15

filtF <- file.path(path, "filtered", paste0(sam, "_F_filt.fastq.gz"))
filtR <- file.path(path, "filtered", paste0(sam, "_R_filt.fastq.gz"))
names(filtF) <- sam
names(filtR) <- sam

FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "GGACTACNVGGGTWTCTAAT"
nchar(c(FWD, REV))

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
primerHits(FWD, fnF[[1]])
primerHits(REV, fnR[[1]])
# Looks like correct primers. Confirmed on visual inspection
# head(substr(unname(getSequences(fnF[[1]])), 1, 30), 20)
# head(substr(unname(getSequences(fnR[[13]])), 1, 30), 20)

out <- filterAndTrim(fnF, filtF, fnR, filtR, 
                     trimLeft=c(19, 20), truncLen=c(200,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
summary(out[,2]/out[,1])
# Median 85% retained
summary(out[,2])
plot(sort(out[,2])) 
# Perhaps 3 failed libraries with very few reads. All others look fine.

plotComplexity(filtF[[13]])
plotComplexity(filtR[[53]])
# Still a low(er) complexity bump. Could be real.

system.time(errF <- learnErrors(filtF, multi=TRUE))
system.time(errR <- learnErrors(filtR, multi=TRUE))
# Like 1 minute each

plotErrors(errF)
plotErrors(errR)
# Looks fine

system.time(ddF <- dada(filtF, err=errF, pool="pseudo", multithread=TRUE)) # ~20 minutes
# saveRDS(ddF, file.path("ddF.rds"))
system.time(ddR <- dada(filtR, err=errR, pool="pseudo", multithread=TRUE)) # ~11 minutes
# saveRDS(ddR, file.path("ddR.rds"))

mm <- mergePairs(ddF, filtF, ddR, filtR, verbose=TRUE)

sta <- makeSequenceTable(mm)
st <- removeBimeraDenovo(sta, method="consensus", multithread=TRUE)
# st <- removeBimeraDenovo(sta, multithread=TRUE) # 
sum(st)/sum(sta)
# 12% chimeras!? Not great, may need to revisit. Was only 6% in Epi1

ncol(sta); ncol(st)
# 90% reduction in ASVs

sq <- getSequences(st)
table(nchar(sq))
# 252-254... like before. Did they change primers since then?
lens <- sort(unique(nchar(sq)))
tot <- sapply(lens, function(l) sum(colSums(st)[nchar(sq)==l]))
names(tot) <- as.character(lens)
tot
# All in 252-254 basically

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(ddF, getN), sapply(ddR, getN), sapply(mm, getN), rowSums(st))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sam
head(track)
# In observed range, but not great, chimera numbers here. Ask Andrea?

tax <- assignTaxonomy(st, "~/tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
saveRDS(st, here("Rdata/following_study","st.rds"))
saveRDS(tax, here("Rdata/following_study","tax.rds"))

system.time(taxsp <- addSpecies(tax, "~/tax/silva_species_assignment_v138.1.fa.gz")) # 8m
saveRDS(taxsp, here("Rdata/following_study","taxsp.rds"))

is.lacto <- taxsp[,6] %in% "Lactobacillus"
sum(is.lacto)
unname(taxsp[is.lacto,6:7])
# Mostly can't speciate lactos this way

# Read in the metadata and coordinate with the sequence table
samdf <- read.csv("epilepsy_mapping.csv")
sapply(samdf, "class")
samdf$Household <- paste0("H", samdf$Household)

# Check if if these match up with the sequence table

head(rownames(st))
head(samdf$SampleID)
head(samdf$SampleID2)

# $SampleID is _J before the sample ID number, and use just the bare number, not the 3-digit version with padded 0s
# $SampleID2 is the character sample Id number.
# This includes controls, and some empty ones, but these aren't in the sequence table...?
samdf[samdf$SampleID2=="",]
samdf[grepl("Cont", samdf$SampleID2),]
# Or maybe they are? Are 94-96 or 122-124 in the SampleID list?
rowSums(st[paste0("Netti0", seq(94,96)),])
# Those are full read-count samples
# 122-124 aren't in the sequenced list
table(samdf$Household, useNA="ifany")
# More to do from here, but let's coordinate first

table(samdf$SampleID)



##foo <- sapply(strsplit(samdf$SampleID, "_"), `[`, 2)
##foo <- gsub("J", "", foo)
##foo <- as.integer(foo)
##foo <- sprintf("%03d", foo)

df <- samdf[!grepl("Cont", samdf$SampleID2),]
df$sam <- sprintf("%03d", as.integer(df$SampleID2))
# Losing info on the 3 controls that were sequences (apparently)
df <- df[df$sam != "0NA",]
dim(df) # 118
df$sam <- paste0("Netti", gsub("Netti", "", df$sam))
# OK...

table(df$sam %in% rownames(st)) # All 118 TRUE
table(rownames(st) %in% df$sam) # 117 TRUE and 1 FALSE !?

rownames(st)[!rownames(st) %in% df$sam] # "Netti051"
# Hm, so one sample not in df, but not the other way?

length(unique(rownames(st))) #118
length(unique(df$sam)) # 117
rownames(st)[!rownames(st) %in% df$sam] # "Netti051"
df$sam[duplicated(df$sam)] # "Netti104"
df[df$sam == "Netti104",] # two entries
# I dunno
# in df, they are both SampleID Netti_J51, and SampleID2 104.

head(unname(taxsp), 20)


