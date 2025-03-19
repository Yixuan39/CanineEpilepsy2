#!/usr/bin/env Rscript

rm(list = ls())
library(tidyverse)
library(tidymodels)
library(phyloseq)
library(here)
library(future)
plan(multisession)

# load data
ps.rare <- readRDS(here('data','following_study','ps_rarefied.rds')) %>% 
    filter_taxa(function(x) sum(x > 0)/length(x) > 0.1, prune = TRUE) %>%
    tax_glom(taxrank = 'Phylum')

sam <- data.frame(sample_data(ps.rare)) %>% 
    dplyr::select(Age..months., Sex, Epileptic.or.Control, Household)
otu <- data.frame(otu_table(ps.rare))

# combine data to a data frame
ps.data <- bind_cols(sam, otu)
performance <- data.frame()

for (i in 1:30) {
    message('Iteration: ', i)
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
    
    # Random Forest model
    rf.model <- rand_forest(mode = "classification",
                            trees = tune(),
                            min_n = tune(),
                            mtry = tune()) %>%
        set_engine("ranger")
    
    
    # Workflow
    rf.workflow <- workflow() %>%
        add_model(rf.model) %>%
        add_recipe(Recipe)
    
    rf.cv <- group_vfold_cv(train_data, group = 'Household', v = 10)
    
    # Tuning grid
    rf.grid <- parameters(trees(),
                          min_n(),
                          mtry = mtry_prop())
    
    # Hyperparameter tuning
    rf.tunning <- rf.workflow %>%
        tune_bayes(
            resamples = rf.cv,
            param_info = rf.grid,
            initial = 5,    
            iter = 40,       
            control = control_bayes(no_improve = 10, save_pred = TRUE)
        )
    
    # Select best hyperparameters
    best.tree <- rf.tunning %>% select_best(metric = "accuracy")
    
    # Fit final model
    rf.fit <- finalize_workflow(rf.workflow, best.tree) %>%
        last_fit(split = data_split)
    
    result <- data.frame(method = 'random forest',
                         data = 'phylum',
                         accuracy = unlist(collect_metrics(rf.fit)[1,3]))
    performance <- performance %>% bind_rows(result)
}

if (!dir.exists(here('data','following_study', 'ml-performance'))) {
    dir.create(here('data','following_study', 'ml-performance'), recursive = TRUE)
}
write_csv(performance, here('data','following_study', 'ml-performance','random-forest-phylum.csv'))

