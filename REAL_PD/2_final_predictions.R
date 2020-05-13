library(plyr)
library(dplyr)
library(glmnet)
library(xgboost)
library(compiler)
enableJIT(3)


#
# Set global parameters
#

# Include ancillary data for training?
incl_ancillary <- FALSE

# Weight training samples in regressions to match evaluation measure?
weighted_regressions <- TRUE


#
# Data should be loaded before running this script
#

# Used features
used_feature_sets <- c('simple_acc_mag_features',
                       'simple_int_acc_mag_features',
                       'simple_gyro_mag_features',
                       'simple_int_gyro_mag_features',
                       #'simple_acc_mag_p_features',
                       #'simple_int_acc_mag_p_features',
                       'ts_features', 
                       'ts_int_features',
                       'ts_gyro_features',
                       'ts_gyro_int_features',
                       #'ts_phone_features',
                       #'ts_phone_int_features',
                       'PMA_features')

# Make data if it does not exist already
if (any(!sapply(used_feature_sets, exists))) {
  stop('You must load data and create features before you can predict!')
}

labels <- merge(labels, type_dat, all = TRUE)

#
# Prepare regression data
#

reg_data <- merge(merge(merge(labels, demographic, all = TRUE), UPDRS3_feat, all = TRUE), UPDRS124, all = TRUE)

for (feat in used_feature_sets) {
  reg_data <- merge(reg_data, get(feat), all = TRUE)
}

reg_data$subject_id <- as.factor(reg_data$subject_id)

# Remove columns with only one level
reg_data <- reg_data[, apply(reg_data, 2, function(x) length(unique(na.omit(x)))) > 1]

# Remove ethnicity
reg_data$Ethnicity <- NULL

# Groups of variables
demo_names <- names(demographic)[-1]
demo_names <- demo_names[demo_names %in% names(reg_data)]

clin_names <- c(names(UPDRS3_feat)[-1], names(UPDRS124)[-1])
clin_names <- clin_names[clin_names %in% names(reg_data)]

for (feat in used_feature_sets) {
  name_var <- paste0(feat, '_names')
  assign(name_var, names(get(feat)))
  assign(name_var, get(name_var)[!grepl('measurement_id', get(name_var))])
  assign(name_var, get(name_var)[get(name_var) %in% names(reg_data)])
}

# Remove ancillary data if specified
if(!incl_ancillary) reg_data <- subset(reg_data, type != 'ancillary')

#
# Results
#

test_predictions <- subset(reg_data, type == 'test')

######################
# On-off best models #
######################

# Prediction data
pred_dat <- subset(reg_data, !is.na(on_off) & !is.na(measurement_id))
test_predictions_on_off <- subset(test_predictions, !(subject_id %in% c('hbv012', 'hbv017', 'hbv018', 'hbv023', 'hbv054')))


# Weights - use local weighting
tmp <- as.data.frame(table(pred_dat$subject_id))
tmp[, 2] <- sqrt(tmp[, 2])
weights <- 1 / tmp[match(pred_dat$subject_id, tmp[, 1]), 2]

#
# Linear models
#

# Null model (everyone not specified below)
on_off_global_subject_mod <- lm(on_off ~ subject_id - 1, data = reg_data[!duplicated(reg_data$measurement_id_original), ])

on_off_global_subject_pred <- predict(on_off_global_subject_mod, newdata = test_predictions_on_off)

# hbv038
f <- formula(paste('on_off ~ subject_id +', paste0(PMA_features_names, collapse = ' + ')))
on_off_PMA_features_subject_lmod <- lm(f, data = pred_dat, weights = weights)

on_off_PMA_features_subject_pred <- predict(on_off_PMA_features_subject_lmod, newdata = test_predictions_on_off)

#
# Subject-specific ridge regression models
#

# hbv043
f <- formula(paste('~ -1 +', paste0(ts_features_names, collapse = ' + ')))
X <- model.matrix(f, pred_dat)
idx_subj <- pred_dat$subject_id == 'hbv043'

on_off_ts_features_per_subject_ridge_mod_hbv043 <- cv.glmnet(x = X[idx_subj,], 
                                                             y = pred_dat$on_off[idx_subj], 
                                                             family = 'gaussian', 
                                                             alpha = 0)

on_off_ts_features_per_subject_ridge_pred_hbv043 <- rep(NA, length = nrow(test_predictions_on_off))
idx_subj <- with(test_predictions_on_off, subject_id == 'hbv043' & !is.na(measurement_id))
on_off_ts_features_per_subject_ridge_pred_hbv043[idx_subj] <- predict(on_off_ts_features_per_subject_ridge_mod_hbv043, 
                                                                      newx = model.matrix(f, test_predictions_on_off[idx_subj, ]),
                                                                      s = 'lambda.min')

# hbv077
f <- formula(paste('~ -1 +', paste0(PMA_features_names, collapse = ' + ')))
X <- model.matrix(f, pred_dat)
idx_subj <- pred_dat$subject_id == 'hbv077'

on_off_PMA_features_per_subject_ridge_mod_hbv077 <- cv.glmnet(x = X[idx_subj,], 
                                                              y = pred_dat$on_off[idx_subj], 
                                                              family = 'gaussian', 
                                                              alpha = 0)

on_off_PMA_features_per_subject_ridge_pred_hbv077 <- rep(NA, length = nrow(test_predictions_on_off))
idx_subj <- with(test_predictions_on_off, subject_id == 'hbv077' & !is.na(measurement_id))
on_off_PMA_features_per_subject_ridge_pred_hbv077[idx_subj] <- predict(on_off_PMA_features_per_subject_ridge_mod_hbv077, 
                                                                       newx = model.matrix(f, test_predictions_on_off[idx_subj, ]),
                                                                       s = 'lambda.min')


# hbv022
f <- formula(paste(' ~ -1 + ', paste0(simple_int_gyro_mag_features_names, collapse = ' + ')))
X <- model.matrix(f, pred_dat)
idx_subj <- pred_dat$subject_id == 'hbv022'

on_off_simple_int_gyro_mag_features_per_subject_ridge_mod_hbv022 <- cv.glmnet(x = X[idx_subj,], 
                                                                              y = pred_dat$on_off[idx_subj], 
                                                                              family = 'gaussian', 
                                                                              alpha = 0)

on_off_simple_int_gyro_mag_features_per_subject_ridge_pred_hbv022 <- rep(NA, length = nrow(test_predictions_on_off))
idx_subj <- with(test_predictions_on_off, subject_id == 'hbv022' & !is.na(measurement_id))
on_off_simple_int_gyro_mag_features_per_subject_ridge_pred_hbv022[idx_subj] <- predict(on_off_simple_int_gyro_mag_features_per_subject_ridge_mod_hbv022, 
                                                                                       newx = model.matrix(f, test_predictions_on_off[idx_subj, ]),
                                                                                       s = 'lambda.min')

# RF model

# hbv014
f <- formula(paste(' ~ subject_id +', paste0(simple_gyro_mag_features_names, collapse = ' + ')))
X <- model.matrix(f, pred_dat)

on_off_simple_gyro_mag_features_subject_rf_xgb_mod <- xgboost(params = list(colsample_bynode = 0.8, 
                                                                            learning_rate = 1, 
                                                                            max_depth = 5, 
                                                                            num_parallel_tree = 500, 
                                                                            subsample = 0.8,
                                                                            nthread = 4), 
                                                              X, 
                                                              weight = weights, 
                                                              label = pred_dat$on_off, 
                                                              nrounds = 1)


on_off_simple_gyro_mag_features_subject_rf_xgb_pred <- rep(NA, length = nrow(test_predictions_on_off))
idx <- with(test_predictions_on_off, !is.na(measurement_id))
on_off_simple_gyro_mag_features_subject_rf_xgb_pred[idx] <- predict(on_off_simple_gyro_mag_features_subject_rf_xgb_mod, 
                                                                    newdata = model.matrix(f, test_predictions_on_off[idx, ]))
#
# Generate test predictions
#

# prediction
on_off <- on_off_global_subject_pred
on_off[test_predictions_on_off$subject_id == 'hbv038'] <- on_off_PMA_features_subject_pred[test_predictions_on_off$subject_id == 'hbv038']
on_off[test_predictions_on_off$subject_id == 'hbv043'] <- on_off_ts_features_per_subject_ridge_pred_hbv043[test_predictions_on_off$subject_id == 'hbv043']
on_off[test_predictions_on_off$subject_id == 'hbv077'] <- on_off_PMA_features_per_subject_ridge_pred_hbv077[test_predictions_on_off$subject_id == 'hbv077']
on_off[test_predictions_on_off$subject_id == 'hbv022'] <- on_off_simple_int_gyro_mag_features_per_subject_ridge_pred_hbv022[test_predictions_on_off$subject_id == 'hbv022']
on_off[test_predictions_on_off$subject_id == 'hbv014'] <- on_off_simple_gyro_mag_features_subject_rf_xgb_pred[test_predictions_on_off$subject_id == 'hbv014']
on_off[is.na(on_off)] <- on_off_global_subject_pred[is.na(on_off)]

on_off[on_off < 0] <- 0
on_off[on_off > 1] <- 1

plot(on_off_global_subject_pred, on_off, asp = 1)

# Fill in predictions
test_predictions$on_off[!(test_predictions$subject_id %in% c('hbv012', 'hbv017', 'hbv018', 'hbv023', 'hbv054'))] <- on_off



##########################
# Dyskinesia best models #
##########################


# Prediction data
pred_dat <- subset(reg_data, !is.na(dyskinesia) & !is.na(measurement_id))
test_predictions_dyskinesia <- subset(test_predictions, !(subject_id %in% c('hbv012', 'hbv014', 'hbv022', 'hbv023', 'hbv038', 'hbv051', 'hbv077')))


# Weights - use local weighting
tmp <- as.data.frame(table(pred_dat$subject_id))
tmp[, 2] <- sqrt(tmp[, 2])
weights <- 1 / tmp[match(pred_dat$subject_id, tmp[, 1]), 2]

#
# Linear models
#

# Null model (everyone not specified below)
dyskinesia_global_subject_mod <- lm(dyskinesia ~ subject_id - 1, data = reg_data[!duplicated(reg_data$measurement_id_original), ])

dyskinesia_global_subject_pred <- predict(dyskinesia_global_subject_mod, newdata = test_predictions_dyskinesia)

#
# Subject-specific ridge regression models
#

# hbv054:
f <- formula(paste('~ -1 + ', paste0(simple_gyro_mag_features_names, collapse = ' + ')))
X <- model.matrix(f, pred_dat)
idx_subj <- pred_dat$subject_id == 'hbv054'

dyskinesia_simple_gyro_mag_features_per_subject_ridge_mod_hbv054 <- cv.glmnet(x = X[idx_subj,], 
                                                                              y = pred_dat$dyskinesia[idx_subj], 
                                                                              family = 'gaussian', 
                                                                              alpha = 0)

dyskinesia_simple_gyro_mag_features_per_subject_ridge_pred_hbv054 <- rep(NA, length = nrow(test_predictions_dyskinesia))
idx_subj <- with(test_predictions_dyskinesia, subject_id == 'hbv054' & !is.na(measurement_id))
dyskinesia_simple_gyro_mag_features_per_subject_ridge_pred_hbv054[idx_subj] <- predict(dyskinesia_simple_gyro_mag_features_per_subject_ridge_mod_hbv054, 
                                                                                       newx = model.matrix(f, test_predictions_dyskinesia[idx_subj, ]),
                                                                                       s = 'lambda.min')


#
# Generate test predictions
#

# prediction
dyskinesia <- dyskinesia_global_subject_pred
dyskinesia[test_predictions_dyskinesia$subject_id == 'hbv054'] <- dyskinesia_simple_gyro_mag_features_per_subject_ridge_pred_hbv054[test_predictions_dyskinesia$subject_id == 'hbv054']
dyskinesia[is.na(dyskinesia)] <- dyskinesia_global_subject_pred[is.na(dyskinesia)]

dyskinesia[dyskinesia < 0] <- 0
dyskinesia[dyskinesia > 2] <- 2

plot(dyskinesia_global_subject_pred, dyskinesia, asp = 1)

# Fill in predictions
test_predictions$dyskinesia[!(test_predictions$subject_id %in% c('hbv012', 'hbv014', 'hbv022', 'hbv023', 'hbv038', 'hbv051', 'hbv077'))] <- dyskinesia


######################
# Tremor best models #
######################


# Prediction data
pred_dat <- subset(reg_data, !is.na(tremor) & !is.na(measurement_id))
test_predictions_tremor <- subset(test_predictions, !(subject_id %in% c('hbv014', 'hbv017', 'hbv018', 'hbv043', 'hbv051', 'hbv077')))


# Weights - use local weighting
tmp <- as.data.frame(table(pred_dat$subject_id))
tmp[, 2] <- sqrt(tmp[, 2])
weights <- 1 / tmp[match(pred_dat$subject_id, tmp[, 1]), 2]

#
# Linear models
#

# Null model (everyone not specified below)
tremor_global_subject_mod <- lm(tremor ~ subject_id - 1, data = reg_data[!duplicated(reg_data$measurement_id_original), ])

tremor_global_subject_pred <- predict(tremor_global_subject_mod, newdata = test_predictions_tremor)



#
# Lasso model
#

# hbv054
f <- formula(paste('~ subject_id + ', paste0(simple_int_acc_mag_features_names, collapse = ' + ')))
X <- model.matrix(f, pred_dat)

tremor_simple_int_acc_mag_features_subject0_lasso_modlambda.min_1  <- cv.glmnet(x = X, 
                                                                                y = pred_dat$tremor,
                                                                                weights = weights,
                                                                                family = 'gaussian', 
                                                                                alpha = 1,
                                                                                relax = TRUE,
                                                                                penalty.factor = as.numeric(!grepl('subject', colnames(X)) | grepl(':', colnames(X))))

tremor_simple_int_acc_mag_features_subject0_lasso_predlambda.min_1 <- rep(NA, length = nrow(test_predictions_tremor))
idx <- with(test_predictions_tremor, !is.na(measurement_id))
tremor_simple_int_acc_mag_features_subject0_lasso_predlambda.min_1[idx] <- predict(tremor_simple_int_acc_mag_features_subject0_lasso_modlambda.min_1, 
                                                                                   newx = model.matrix(f, test_predictions_tremor[idx, ]),
                                                                                   s = 'lambda.min',
                                                                                   gamma = 1)

#
# Subject-specific ridge/lasso regression models
#

# hbv038:
f <- formula(paste('~ ', paste0(ts_features_names, collapse = ' + ')))
X <- model.matrix(f, pred_dat)
idx_subj <- pred_dat$subject_id == 'hbv038'

tremor_ts_features_per_subject_ridge_mod_hbv038 <- cv.glmnet(x = X[idx_subj,], 
                                                             y = pred_dat$tremor[idx_subj], 
                                                             family = 'gaussian', 
                                                             alpha = 0)

tremor_ts_features_per_subject_ridge_pred_hbv038 <- rep(NA, length = nrow(test_predictions_tremor))
idx_subj <- with(test_predictions_tremor, subject_id == 'hbv038' & !is.na(measurement_id))
tremor_ts_features_per_subject_ridge_pred_hbv038[idx_subj] <- predict(tremor_ts_features_per_subject_ridge_mod_hbv038, 
                                                                      newx = model.matrix(f, test_predictions_tremor[idx_subj, ]),
                                                                      s = 'lambda.min')


# hbv022:
f <- formula(paste('~ ', paste0(PMA_features_names, collapse = ' + ')))
X <- model.matrix(f, pred_dat)
idx_subj <- pred_dat$subject_id == 'hbv022'

tremor_PMA_features_per_subject_ridge_mod_hbv022 <- cv.glmnet(x = X[idx_subj,], 
                                                              y = pred_dat$tremor[idx_subj], 
                                                              family = 'gaussian', 
                                                              alpha = 0)

tremor_PMA_features_per_subject_ridge_pred_hbv022 <- rep(NA, length = nrow(test_predictions_tremor))
idx_subj <- with(test_predictions_tremor, subject_id == 'hbv022' & !is.na(measurement_id))
tremor_PMA_features_per_subject_ridge_pred_hbv022[idx_subj] <- predict(tremor_PMA_features_per_subject_ridge_mod_hbv022 , 
                                                                       newx = model.matrix(f, test_predictions_tremor[idx_subj, ]),
                                                                       s = 'lambda.min')


# hbv023:
f <- formula(paste('~ -1 + ', paste0(unlist(sapply(paste0(used_feature_sets, '_names'), get)), collapse = ' + ')))
X <- model.matrix(f, pred_dat)
idx_subj <- pred_dat$subject_id == 'hbv023'

tremor_all_per_subject_lasso_modlambda.1se_1_hbv023 <- cv.glmnet(x = X[idx_subj,], 
                                                                 y = pred_dat$tremor[idx_subj], 
                                                                 family = 'gaussian', 
                                                                 alpha = 1,
                                                                 relax = TRUE,
                                                                 pmax = 200)

tremor_all_per_subject_lasso_predlambda.1se_1_hbv023 <- rep(NA, length = nrow(test_predictions_tremor))
idx_subj <- with(test_predictions_tremor, subject_id == 'hbv023' & !is.na(measurement_id))
tremor_all_per_subject_lasso_predlambda.1se_1_hbv023[idx_subj] <- predict(tremor_all_per_subject_lasso_modlambda.1se_1_hbv023 , 
                                                                          newx = model.matrix(f, test_predictions_tremor[idx_subj, ]),
                                                                          s = 'lambda.1se',
                                                                          gamma = 1)



# RF model

# hbv012
f <- formula(paste(' ~ subject_id +', paste0(simple_gyro_mag_features_names, collapse = ' + ')))
X <- model.matrix(f, pred_dat)

tremor_simple_int_gyro_mag_features_subject_rf_xgb_mod  <- xgboost(params = list(colsample_bynode = 0.8, 
                                                                                 learning_rate = 1, 
                                                                                 max_depth = 5, 
                                                                                 num_parallel_tree = 500, 
                                                                                 subsample = 0.8,
                                                                                 nthread = 4), 
                                                                   X, 
                                                                   weight = weights, 
                                                                   label = pred_dat$tremor, 
                                                                   nrounds = 1)


tremor_simple_int_gyro_mag_features_subject_rf_xgb_pred <- rep(NA, length = nrow(test_predictions_tremor))
idx <- with(test_predictions_tremor, !is.na(measurement_id))
tremor_simple_int_gyro_mag_features_subject_rf_xgb_pred[idx] <- predict(tremor_simple_int_gyro_mag_features_subject_rf_xgb_mod , 
                                                                         newdata = model.matrix(f, test_predictions_tremor[idx, ]))


#
# Generate test predictions
#

# prediction
tremor <- tremor_global_subject_pred
tremor[test_predictions_tremor$subject_id == 'hbv054'] <- tremor_simple_int_acc_mag_features_subject0_lasso_predlambda.min_1[test_predictions_tremor$subject_id == 'hbv054']
tremor[test_predictions_tremor$subject_id == 'hbv038'] <- tremor_ts_features_per_subject_ridge_pred_hbv038[test_predictions_tremor$subject_id == 'hbv038']
tremor[test_predictions_tremor$subject_id == 'hbv022'] <- tremor_PMA_features_per_subject_ridge_pred_hbv022[test_predictions_tremor$subject_id == 'hbv022']
tremor[test_predictions_tremor$subject_id == 'hbv023'] <- tremor_all_per_subject_lasso_predlambda.1se_1_hbv023[test_predictions_tremor$subject_id == 'hbv023']
tremor[test_predictions_tremor$subject_id == 'hbv012'] <- tremor_simple_int_gyro_mag_features_subject_rf_xgb_pred[test_predictions_tremor$subject_id == 'hbv012']
tremor[is.na(tremor)] <- tremor_global_subject_pred[is.na(tremor)]

tremor[tremor < 0] <- 0
tremor[tremor > 4] <- 4

plot(tremor_global_subject_pred, tremor, asp = 1)

# Fill in predictions
test_predictions$tremor[!(test_predictions$subject_id %in% c('hbv014', 'hbv017', 'hbv018', 'hbv043', 'hbv051', 'hbv077'))] <- tremor



#
# Test prediction
#

# Aggregate predictions for individuals with several devices
test_predictions <- test_predictions[, c('measurement_id_original', 'on_off', 'dyskinesia', 'tremor')]
names(test_predictions)[1] <- 'measurement_id'
test_predictions <- aggregate(. ~ measurement_id, test_predictions, FUN = mean, na.rm = TRUE, na.action = NULL)

#
# Test prediction
#

OnOff_Submission <- test_predictions[, c('measurement_id', 'on_off')]
names(OnOff_Submission)[2] <- 'prediction'

Dyskinesia_Submission <- test_predictions[, c('measurement_id', 'dyskinesia')]
names(Dyskinesia_Submission)[2] <- 'prediction'

Tremor_Submission <- test_predictions[, c('measurement_id', 'tremor')]
names(Tremor_Submission)[2] <- 'prediction'

