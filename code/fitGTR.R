library(tidyverse)
library(phyloseq)
library(vegan)
library(phangorn)
library(DECIPHER)
library(here)

filt_fun <- function (x, min_reads = 100, min_samples = 5) {
    (sum(x) > min_reads) & (sum(x > 0) > min_samples)
}

ps <- readRDS(here('Rdata','following_study','ps.rds')) %>% 
    filter_taxa(filt_fun, prune = TRUE) %>%
    transform_sample_counts(function(x) x / sum(x) )

alignment <- AlignSeqs(refseq(ps), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

saveRDS(fitGTR, here('Rdata','following_study','fitGTR_filtered.rds'))