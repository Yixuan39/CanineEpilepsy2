library(tidyverse)
library(phyloseq)
library(vegan)
library(phangorn)
library(DECIPHER)
library(here)

ps <- readRDS(here('Rdata','following_study','ps.rds')) %>%
    transform_sample_counts(function(x) x / sum(x) )

alignment <- AlignSeqs(refseq(ps), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

saveRDS(fitGTR, here('Rdata','following_study','fitGTR.rds'))