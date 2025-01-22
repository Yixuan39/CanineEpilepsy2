#!/usr/bin/env Rscript
rm(list = ls())
set.seed(2024)
library(tidyverse)
library(tidymodels)
library(here)
library(vip)
library(future)
plan(multisession)
# load data
ps.data <- readRDS(here('data','following_study','ml_data.rds'))

# Split data
data_split <- initial_split(ps.data, prop = 0.7, strata = Epileptic.or.Control)
train_data <- training(data_split)
test_data <- testing(data_split)

# Recipe
Recipe <- recipe(Epileptic.or.Control ~ ., data = train_data) %>%
    step_zv(all_predictors()) %>%
    step_corr(all_numeric_predictors(), threshold = 0.9) %>%
    step_dummy(Sex)

# Model 
log.model <- logistic_reg(penalty = tune(), mixture = tune()) %>%
    set_mode('classification') %>% 
    set_engine("glmnet")

# Workflow
log.workflow <- workflow() %>%
    add_model(log.model) %>%
    add_recipe(Recipe)

# Grid
log.grid <- grid_random(penalty(), mixture(), size = 1000)

# Hyperparameter tuning
log.tunning <- log.workflow %>%
    tune_grid(
        resamples = bootstraps(train_data, times = 1000),
        grid = log.grid,
        metrics = metric_set(roc_auc, pr_auc, accuracy),
        control = control_grid(save_pred = TRUE)
    )

best.log <- log.tunning %>% select_best(metric = "accuracy")

# Finalize workflow
final.log.workflow <- finalize_workflow(log.workflow, best.log)

log.fit <- final.log.workflow %>%
    last_fit(split = data_split)

log.fit %>% collect_metrics()

# Random Forest model
rf.model <- rand_forest(mode = "classification",
                        trees = tune(),
                        min_n = tune(),
                        mtry = tune()) %>%
    set_engine("ranger", importance = "impurity") %>%
    set_mode("classification")


# Workflow
rf.workflow <- workflow() %>%
    add_model(rf.model) %>%
    add_recipe(Recipe)

rf.boot <- bootstraps(train_data, times = 1000)

# Tuning grid
rf.grid <- grid_random(
    trees(),min_n(),
    finalize(mtry(), rf.boot),
    size = 1000
)

# Hyperparameter tuning
rf.tunning <- rf.workflow %>%
    tune_grid(
        resamples = rf.boot,
        grid = rf.grid,
        metrics = metric_set(roc_auc, pr_auc, accuracy),
        control = control_grid(save_pred = TRUE)
    )

# Select best hyperparameters
best.tree <- rf.tunning %>% select_best(metric = "accuracy")

# Finalize workflow
final.rf.workflow <- finalize_workflow(rf.workflow, best.tree)

# Fit final model
rf.fit <- final.rf.workflow %>%
    last_fit(split = data_split)

rf.fit %>%
    collect_metrics()

tree <- extract_workflow(rf.fit)

tree %>% 
    extract_fit_parsnip() %>% 
    vip()

# Random Forest model
svm.model <- svm_linear(cost = tune()) %>% 
    set_engine("kernlab") %>%
    set_mode("classification")


# Workflow
svm.workflow <- workflow() %>%
    add_model(svm.model) %>%
    add_recipe(Recipe)

# Tuning grid
svm.grid <- grid_random(cost(), size = 1000)

# Hyperparameter tuning
svm.tunning <- svm.workflow %>%
    tune_grid(
        resamples = bootstraps(train_data, times = 1000),
        grid = svm.grid,
        metrics = metric_set(roc_auc, pr_auc, accuracy),
        control = control_grid(save_pred = TRUE)
    )

# Select best hyperparameters
best.svm <- svm.tunning %>% select_best(metric = "accuracy")

# Finalize workflow
final.svm.workflow <- finalize_workflow(svm.workflow, best.svm)

# Fit final model
svm.fit <- final.svm.workflow %>%
    last_fit(split = data_split)

svm.fit %>%
    collect_metrics()

# Random Forest model
mlp.model <- mlp(hidden_units = tune(),
                 penalty = tune(),
                 epochs = tune()) %>% 
    set_engine("nnet") %>%
    set_mode("classification") %>% 
    set_args(MaxNWts = 100000)


# Workflow
mlp.workflow <- workflow() %>%
    add_model(mlp.model) %>%
    add_recipe(Recipe)

# Tuning grid
mlp.grid <- grid_random(hidden_units(), penalty(), epochs(), size = 1000)

# Hyperparameter tuning
mlp.tunning <- mlp.workflow %>%
    tune_grid(
        resamples = bootstraps(train_data, times = 1000),
        grid = mlp.grid,
        metrics = metric_set(roc_auc, pr_auc, accuracy),
        control = control_grid(save_pred = TRUE)
    )

# Select best hyperparameters
best.mlp <- mlp.tunning %>% select_best(metric = "accuracy")

# Finalize workflow
final.mlp.workflow <- finalize_workflow(mlp.workflow, best.mlp)

# Fit final model
mlp.fit <- final.mlp.workflow %>%
    last_fit(split = data_split)

mlp.fit %>%
    collect_metrics()

# save environment
save.image(here('data','following_study','ml_models.rds'))
