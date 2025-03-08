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

pcoa <- ps.rare %>% 
    ordinate(method = 'PCoA', distance = 'bray')
n.vecter <- 30

sam <- data.frame(sample_data(ps.rare)) %>% 
    dplyr::select(Age..months., Sex, Epileptic.or.Control, Household)

# combine data to a data frame
ps.data <- bind_cols(sam, pcoa$vectors[,1:n.vecter])
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
    
    # set up model
    svm.model <- svm_rbf(cost = tune(), rbf_sigma = tune()) %>% 
        set_engine("kernlab") %>%
        set_mode("classification")
    
    
    # Workflow
    svm.workflow <- workflow() %>%
        add_model(svm.model) %>%
        add_recipe(Recipe)
    
    # Tuning grid
    svm.grid <- grid_regular(cost(),
                             rbf_sigma(),
                             levels = c(cost = 10, rbf_sigma = 10))
    
    # Hyperparameter tuning
    svm.tunning <- svm.workflow %>%
        tune_grid(
            resamples = group_vfold_cv(train_data, group = 'Household', v = 10, repeats = 10),
            grid = svm.grid,
            control = control_grid(save_pred = TRUE)
        )
    
    # Select best hyperparameters
    best.svm <- svm.tunning %>% select_best(metric = "accuracy")
    
    svm.fit <- finalize_workflow(svm.workflow, best.svm) %>%
        last_fit(split = data_split)
    
    result <- data.frame(method = 'SVM',
                         data = 'PCoA',
                         accuracy = unlist(collect_metrics(svm.fit)[1,3]))
    performance <- performance %>% bind_rows(result)
}

if (!dir.exists(here('data','following_study', 'ml-performance'))) {
    dir.create(here('data','following_study', 'ml-performance'), recursive = TRUE)
}
write_csv(performance, here('data','following_study', 'ml-performance','support-vector-machine-pcoa.csv'))

