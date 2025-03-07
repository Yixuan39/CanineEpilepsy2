#!/usr/bin/env Rscript

rm(list = ls())
library(tidyverse)
library(tidymodels)
library(phyloseq)
library(here)
library(future)
plan(multisession)
theme_set(theme_bw())
# load data
ps.rare <- readRDS(here('data','following_study','ps_rarefied.rds')) %>% 
    filter_taxa(function(x) sum(x > 0)/length(x) > 0.1, prune = TRUE)

sam <- data.frame(sample_data(ps.rare)) %>% 
    dplyr::select(Age..months., Sex, Epileptic.or.Control, Household)
otu <- data.frame(otu_table(ps.rare))

# combine data to a data frame
ps.data <- bind_cols(sam, otu)
performance <- data.frame()

for (i in 1:30) {
    set.seed(i)
    # Split data
    data_split <- group_initial_split(ps.data, group = 'Household', prop = 0.7)
    train_data <- training(data_split)
    test_data <- testing(data_split)
    
    # Recipe
    Recipe <- recipe(Epileptic.or.Control ~ ., data = train_data) %>%
        update_role(Household, new_role = 'group variable') %>%
        step_nzv(all_predictors()) %>%
        step_dummy(all_nominal_predictors()) %>%
        step_normalize(Age..months., starts_with('ASV'))
    
    # Model 
    log.model <- logistic_reg(penalty = tune(), mixture = tune()) %>%
        set_mode('classification') %>% 
        set_engine("glmnet")
    
    # Workflow
    log.workflow <- workflow() %>%
        add_model(log.model) %>%
        add_recipe(Recipe)
    
    # Grid
    log.grid <- grid_regular(penalty(), mixture(),
                             levels = c(penalty = 10, mixture = 10))
    
    # Hyperparameter tuning
    log.tunning <- log.workflow %>%
        tune_grid(
            resamples = group_vfold_cv(train_data, group = 'Household', v = 10, repeats = 10),
            grid = log.grid,
            control = control_grid(save_pred = TRUE)
        )
    
    best.log <- log.tunning %>% select_best(metric = "accuracy")
    
    log.fit <- finalize_workflow(log.workflow, best.log) %>%
        last_fit(split = data_split)
    
    result <- data.frame(method = 'elastic net',
                         data = 'ASV',
                         accuracy = unlist(collect_metrics(log.fit)[1,3]))
    performance <- performance %>% bind_rows(result)
}

if (!dir.exists(here('data','following_study', 'ml-performance'))) {
    dir.create(here('data','following_study', 'ml-performance'), recursive = TRUE)
}
write_csv(performance, here('data','following_study', 'ml-performance','logistic-ASV.csv'))

