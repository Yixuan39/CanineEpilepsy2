#!/usr/bin/env Rscript
set.seed(2024)
library(tidyverse)
library(tidymodels)
library(here)
library(future)
plan(multisession)
# load data
data_split <- readRDS(here('data','following_study','data_split.rds'))

# Split data
train_data <- training(data_split)
test_data <- testing(data_split)

# Recipe
Recipe <- recipe(Epileptic.or.Control ~ ., data = train_data) %>%
    step_zv(all_predictors()) %>%
    step_corr(all_numeric_predictors(), threshold = 0.9) %>%
    step_dummy(Sex)

#####################
# Elastic Net model #
#####################

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

# tuning
log.tunning <- log.workflow %>%
    tune_grid(
        resamples = bootstraps(train_data, times = 1000),
        grid = log.grid,
        metrics = metric_set(roc_auc, pr_auc, accuracy),
        control = control_grid(save_pred = TRUE)
    )

best.log <- log.tunning %>% select_best(metric = "accuracy")

log.fit <- finalize_workflow(log.workflow, best.log) %>%
    last_fit(split = data_split)

log.fit %>% collect_metrics()


#################
# Random Forest #
#################

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

# Grid
rf.grid <- grid_random(
    trees(),min_n(),
    finalize(mtry(), rf.boot),
    size = 1000
)

# tuning
rf.tunning <- rf.workflow %>%
    tune_grid(
        resamples = rf.boot,
        grid = rf.grid,
        metrics = metric_set(roc_auc, pr_auc, accuracy),
        control = control_grid(save_pred = TRUE)
    )

# Select best hyperparameters
best.tree <- rf.tunning %>% select_best(metric = "accuracy")

# Fit final model
rf.fit <- finalize_workflow(rf.workflow, best.tree) %>%
    last_fit(split = data_split)

rf.fit %>% collect_metrics()

##########################
# Support Vector Machine #
##########################

svm.model <- svm_rbf(cost = tune(), rbf_sigma = tune()) %>% 
    set_engine("kernlab") %>%
    set_mode("classification")


# Workflow
svm.workflow <- workflow() %>%
    add_model(svm.model) %>%
    add_recipe(Recipe)

# Tuning grid
svm.grid <- grid_random(cost(), rbf_sigma(), size = 1000)

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

# Fit final model
svm.fit <- finalize_workflow(svm.workflow, best.svm) %>%
    last_fit(split = data_split)

svm.fit %>% collect_metrics()


#########################
# Multilayer Perceptron #
#########################

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

# Fit final model
mlp.fit <- finalize_workflow(mlp.workflow, best.mlp) %>%
    last_fit(split = data_split)

mlp.fit %>% collect_metrics()

session.info <- sessioninfo::session_info()

final.models <- list(elastic.net = log.fit,
                     random.forest = rf.fit,
                     svm = svm.fit, 
                     mlp = mlp.fit,
                     session.info = session.info)

# save environment
saveRDS(here('data','following_study','ml_models.rds'))
