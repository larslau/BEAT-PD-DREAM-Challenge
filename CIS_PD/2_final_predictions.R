library(plyr)
library(dplyr)
library(glmnet)
library(xgboost)
library(pliable)

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
                       'ts_features', 
                       'ts_int_features',
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
  assign(name_var, names(get(feat))[-1])
  assign(name_var, get(name_var)[get(name_var) %in% names(reg_data)])
}

# Remove ancillary data if specified
if(!incl_ancillary) reg_data <- subset(reg_data, type != 'ancillary')


#
# Custom CV function for pliable lasso
#

cv.pliable <- function (fit, x, z, y, nfolds = 10, foldid = NULL, keep = F, 
                        type.measure = c("deviance", "class"), verbose = TRUE) {
  status = NULL
  if (fit$family == "cox") {
    o = order(y - 1e-04 * status)
    x = x[o, ]
    z = z[o, ]
    y = y[o]
    status = status[o]
  }
  type.measure = match.arg(type.measure)
  if (fit$family == "gaussian") 
    errfun = pliable:::errfun.gaussian
  if (fit$family == "binomial" & type.measure == "deviance") 
    errfun = pliable:::errfun.binomial
  if (fit$family == "binomial" & type.measure == "class") 
    errfun = function(y, yhat, w) 1 * (y != yhat) * w
  if (fit$family == "cox" & is.null(status)) 
    stop("Must supply status arg with family='cox'")
  BIG = 1e+10
  ni = ncol(x)
  no = length(y)
  nz = ncol(z)
  if (fit$family == "cox") 
    coxinfo.folds = vector("list", nfolds)
  ggg = vector("list", nfolds)
  yhat = array(NA, c(no, length(fit$lambda)))
  if (is.null(foldid)) 
    foldid = sample(rep(seq(nfolds), length = no))
  nfolds = length(table(foldid))
  status.in = NULL
  for (ii in 1:nfolds) {
    oo = foldid == ii
    if (fit$family == "cox") 
      status.in = status[!oo]
    if (verbose) 
      cat(c("\n", "FOLD=", ii), fill = T)
    ggg[[ii]] = pliable(x[!oo, , drop = F], z[!oo, , drop = F], 
                        y[!oo], family = fit$family, lambda = fit$lambda, 
                        alpha = fit$args$alpha, lambda.min.ratio = fit$args$lambda.min.ratio, 
                        w = fit$w[!oo], thr = fit$args$thr, maxit = fit$args$maxit, 
                        tt = fit$args$tt, mxthit = fit$args$mxthit, mxkbt = fit$args$mxkbt, 
                        mxth = fit$args$mxth, kpmax = fit$args$kpmax, kthmax = fit$args$kthmax, 
                        maxinter = 10000,
                        screen = fit$args$screen)
    if (fit$family == "gaussian" | fit$family == "binomial") 
      yhat[oo, ] = predict.pliable(ggg[[ii]], x[oo, , drop = F], 
                                   z[oo, ])
    if (fit$family == "cox") 
      coxinfo.folds[[ii]] = ggg[[ii]]$coxinfo
  }
  nonzero = colSums(fit$beta != 0)
  yhat.preval = NULL
  ym = array(y, dim(yhat))
  if (type.measure == "class") 
    yhat = 1 * (yhat > mean(y))
  if (fit$family == "gaussian" | fit$family == "binomial") {
    err = errfun(ym, yhat, fit$w)
    cvm = apply(err, 2, mean, na.rm = T)
    nn = apply(!is.na(err), 2, sum, na.rm = T)
    cvsd = sqrt(apply(err, 2, var, na.rm = T)/nn)
    if (keep) 
      yhat.preval = yhat
  }
  if (fit$family == "cox") {
    err = matrix(NA, nfolds, length(fit$lambda))
    for (ii in 1:nfolds) {
      oo = foldid == ii
      fit1 = predict.pliable(ggg[[ii]], x, z)
      fit2 = predict.pliable(ggg[[ii]], x[!oo, , drop = F], 
                             z[!oo, ])
      for (k in 1:length(fit$lambda)) {
        dev1 = devc(no, fit$coxinfo$kq, fit$coxinfo$ddq, 
                    fit1[, k], fit$coxinfo$iriskq, status)
        dev2 = devc(sum(!oo), coxinfo.folds[[ii]]$kq, 
                    coxinfo.folds[[ii]]$ddq, fit2[, k], coxinfo.folds[[ii]]$iriskq, 
                    status[!oo])
        err[ii, k] = dev1 - dev2
      }
    }
    cvm = apply(err, 2, mean, na.rm = T)
    nn = apply(!is.na(err), 2, sum, na.rm = T)
    cvsd = sqrt(apply(err, 2, var, na.rm = T)/nn)
  }
  cvm.nz = cvm
  cvm.nz[nonzero == 0] = BIG
  imin = which.min(cvm.nz)
  imin.1se = which(cvm < cvm[imin] + cvsd[imin])[1]
  out = list(lambda = fit$lambda, cvm = cvm, cvsd = cvsd, cvup = cvm + 
               cvsd, cvlo = cvm - cvsd, nz = nonzero, df = nonzero, 
             yhat.preval = yhat.preval, lambda.min = fit$lambda[imin], 
             lambda.1se = fit$lambda[imin.1se], name = "Error")
  class(out) = "cv.pliable"
  return(out)
}




#
# Results
#

test_predictions <- subset(reg_data, type == 'test')

######################
# On-off best models #
######################


# "on_off_ts_features_subject0_pliable_0.5_mod" - 1038  
# 
# "on_off_all_per_subject_lasso_modlambda.min_gamma.min" - Seems unstable
# 
# "on_off_ts_int_features_subject0_lasso_modlambda.min_1" - 1044 
#
# "on_off_simple_int_acc_mag_features_subject0_lasso_modlambda.min_1" - 1032 
#
# "on_off_comb_per_subject_lasso_modlambda.min_1" - 1006 
#
# "on_off_all_subject0_pliable_0.5_mod" - 1020, 1051, 1048, 1049   
#
# "on_off_ts_features_subject0_ridge_mod" - 1007  
#
# "on_off_simple_acc_mag_features_subject0_ridge_mod" - 1039  
#
# "on_off_comb_subject0_ridge_mod" - 1004  
#
# "on_off_ts_int_features_subject0_lasso_modlambda.min_gamma.min" - 1019, 1043    
# 
# "on_off_comb_subject0_lasso_modlambda.min_gamma.min" - 1034

# Training data
idx_train <- with(reg_data, !is.na(on_off) & type != 'test')
train_dat <- subset(reg_data, idx_train)

# Weights
if (weighted_regressions) {
  tmp <- as.data.frame(table(reg_data$subject_id))
  tmp[, 2] <- sqrt(tmp[, 2])
  weights <- 1 / tmp[match(subset(reg_data, !is.na(on_off))$subject_id, tmp[, 1]), 2]
} else {
  weights <- NULL
}

# 1. model: null model
mod1 <- lm(on_off ~ subject_id, data = train_dat, weights = weights)
pred1 <- predict(mod1, newdata = subset(test_predictions, subject_id %in% train_dat$subject_id))

# 2. model
# on_off_ts_features_subject0_pliable_0.5_mod
f <- formula(paste('~ subject_id +', paste0(ts_features_names, collapse = ' + ')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

idx_keep <- apply(X, 2, sd) != 0
X <- X[, idx_keep]
X_test <- X_test[, idx_keep]

non_subj_cols <- !grepl('subject', colnames(X)) | grepl(':', colnames(X))

mx <- colMeans(X[, non_subj_cols])
sx <- sqrt(apply(X[, non_subj_cols], 2, var))
X[, non_subj_cols] <- scale(X[, non_subj_cols], mx, sx) 
X_test[, non_subj_cols] <- scale(X_test[, non_subj_cols], mx, sx) 


# Z design matrix
Z <- model.matrix(~ subject_id, train_dat)
Z_test <- model.matrix(~ subject_id, subset(test_predictions, subject_id %in% train_dat$subject_id))

idx_keep <- apply(Z, 2, sd) != 0
Z <- Z[, idx_keep]
Z_test <- Z_test[, idx_keep]

# Penalty factor
penalty.factor <- as.numeric(!grepl('subject', colnames(X)) | grepl(':', colnames(X)))

# Pliable lasso model
mod_pl <- pliable(x = X, 
                  z = Z, 
                  y = train_dat$on_off, 
                  alpha = 0.5,
                  w = weights,
                  family = 'gaussian', 
                  penalty.factor = penalty.factor, 
                  zlinear = FALSE,
                  maxinter = 2000)

mod_pl_cv <- cv.pliable(mod_pl, 
                        x = X, 
                        z = Z, 
                        y = train_dat$on_off, 
                        nfolds = 5,
                        verbose = FALSE)

pred2 <- predict(mod_pl, 
                 x = X_test, 
                 z = Z_test, 
                 lambda = mod_pl_cv$lambda.min)



# 3. model
# on_off_all_per_subject_lasso_modlambda.min_gamma.min

f <- formula(paste('~ -1 + ', paste0(unlist(sapply(paste0(used_feature_sets, '_names'), get)), collapse = ' + ')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

mod3 <- list()
pred3 <- rep(NA, nrow(X_test))

for (subj in unique(train_dat$subject_id)) {
  # Fit model and predict on subject data
  idx_subj <- with(train_dat, subject_id == subj)
  mod3[[subj]] <- cv.glmnet(x = X[idx_subj,], y = train_dat$on_off[idx_subj], 
                            family = 'gaussian', 
                            alpha = 1, 
                            pmax = 200, 
                            relax = TRUE)
  idx_subj <- with(subset(test_predictions, subject_id %in% train_dat$subject_id), subject_id == subj)
  
  pred3[idx_subj] <- predict(mod3[[subj]], newx = X_test[idx_subj, ], s = 'lambda.min', gamma = 'gamma.min')
}

# 4. model
# on_off_ts_int_features_subject0_lasso_modlambda.min_1
f <- formula(paste('~ subject_id +', paste0(ts_int_features_names, collapse = ' + ')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

mod4 <- cv.glmnet(x = X, 
                  y = train_dat$on_off, 
                  weights = weights, 
                  family = 'gaussian',
                  penalty.factor = as.numeric(!grepl('subject', colnames(X)) | grepl(':', colnames(X))), 
                  alpha = 1, 
                  pmax = 200, 
                  relax = TRUE)

pred4 <- predict(mod4, newx = X_test, s = 'lambda.min', gamma = 1)

# 5. model
# on_off_simple_int_acc_mag_features_subject0_lasso_modlambda.min_1
f <- formula(paste('~ subject_id +', paste0(simple_int_acc_mag_features_names, collapse = ' + ')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

mod5 <- cv.glmnet(x = X, 
                  y = train_dat$on_off, 
                  weights = weights, 
                  family = 'gaussian',
                  penalty.factor = as.numeric(!grepl('subject', colnames(X)) | grepl(':', colnames(X))), 
                  alpha = 1, 
                  pmax = 200, 
                  relax = TRUE)

pred5 <- predict(mod5, newx = X_test, s = 'lambda.min', gamma = 1)

# 6. model
# on_off_comb_per_subject_lasso_modlambda.min_1

f <- formula(paste('~ -1 +', paste0(simple_int_acc_mag_features_names, collapse = ' + '), '+', 
                   paste0(ts_features_names, collapse = ' + ')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

mod6 <- list()
pred6 <- rep(NA, nrow(X_test))

for (subj in unique(train_dat$subject_id)) {
  # Fit model and predict on subject data
  idx_subj <- with(train_dat, subject_id == subj)
  mod6[[subj]] <- cv.glmnet(x = X[idx_subj,], y = train_dat$on_off[idx_subj], 
                            family = 'gaussian', 
                            alpha = 1, 
                            pmax = 200, 
                            relax = TRUE)
  idx_subj <- with(subset(test_predictions, subject_id %in% train_dat$subject_id), subject_id == subj)
  
  pred6[idx_subj] <- predict(mod6[[subj]], newx = X_test[idx_subj, ], s = 'lambda.min', gamma = 1)
}


# 7. model
# on_off_all_subject0_pliable_0.5_mod
f <- formula(paste('~ subject_id +', paste0(unlist(sapply(paste0(used_feature_sets, '_names'), get)), collapse = ' + ')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

idx_keep <- apply(X, 2, sd) != 0
X <- X[, idx_keep]
X_test <- X_test[, idx_keep]

non_subj_cols <- !grepl('subject', colnames(X)) | grepl(':', colnames(X))

mx <- colMeans(X[, non_subj_cols])
sx <- sqrt(apply(X[, non_subj_cols], 2, var))
X[, non_subj_cols] <- scale(X[, non_subj_cols], mx, sx) 
X_test[, non_subj_cols] <- scale(X_test[, non_subj_cols], mx, sx) 


# Z design matrix
Z <- model.matrix(~ subject_id, train_dat)
Z_test <- model.matrix(~ subject_id, subset(test_predictions, subject_id %in% train_dat$subject_id))

idx_keep <- apply(Z, 2, sd) != 0
Z <- Z[, idx_keep]
Z_test <- Z_test[, idx_keep]

# Penalty factor
penalty.factor <- as.numeric(!grepl('subject', colnames(X)) | grepl(':', colnames(X)))

# Pliable lasso model
mod_pl <- pliable(x = X, 
                  z = Z, 
                  y = train_dat$on_off, 
                  alpha = 0.5,
                  w = weights,
                  family = 'gaussian', 
                  penalty.factor = penalty.factor, 
                  zlinear = FALSE,
                  maxinter = 2000)

mod_pl_cv <- cv.pliable(mod_pl, 
                        x = X, 
                        z = Z, 
                        y = train_dat$on_off, 
                        nfolds = 5,
                        verbose = FALSE)

pred7 <- predict(mod_pl, 
                 x = X_test, 
                 z = Z_test, 
                 lambda = mod_pl_cv$lambda.min)


# 8. model
# on_off_ts_features_subject0_ridge_mod
f <- formula(paste('~ subject_id +', paste0(ts_features_names, collapse = ' + ')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

mod8 <- cv.glmnet(x = X, 
                  y = train_dat$on_off, 
                  weights = weights, 
                  family = 'gaussian',
                  penalty.factor = as.numeric(!grepl('subject', colnames(X)) | grepl(':', colnames(X))), 
                  alpha = 0)

pred8 <- predict(mod8, newx = X_test, s = 'lambda.min')

# 9. model
# on_off_simple_acc_mag_features_subject0_ridge_mod
f <- formula(paste('~ subject_id +', paste0(simple_acc_mag_features_names, collapse = ' + ')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

mod9 <- cv.glmnet(x = X, 
                  y = train_dat$on_off, 
                  weights = weights, 
                  family = 'gaussian',
                  penalty.factor = as.numeric(!grepl('subject', colnames(X)) | grepl(':', colnames(X))), 
                  alpha = 0)

pred9 <- predict(mod9, newx = X_test, s = 'lambda.min')

# 10. model
# on_off_comb_subject0_ridge_mod
f <- formula(paste('~ subject_id +', paste0(simple_int_acc_mag_features_names, collapse = ' + '), '+',
                   paste0(ts_features_names, collapse = ' + ')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

mod10 <- cv.glmnet(x = X, 
                   y = train_dat$on_off, 
                   weights = weights, 
                   family = 'gaussian',
                   penalty.factor = as.numeric(!grepl('subject', colnames(X)) | grepl(':', colnames(X))), 
                   alpha = 0)

pred10 <- predict(mod10, newx = X_test, s = 'lambda.min')

# 11. model
# on_off_ts_int_features_subject0_lasso_modlambda.min_gamma.min
f <- formula(paste('~ subject_id +', paste0(ts_int_features_names, collapse = ' + ')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

mod11 <- cv.glmnet(x = X, 
                   y = train_dat$on_off, 
                   weights = weights, 
                   family = 'gaussian',
                   penalty.factor = as.numeric(!grepl('subject', colnames(X)) | grepl(':', colnames(X))), 
                   alpha = 1, 
                   pmax = 200, 
                   relax = TRUE)

pred11 <- predict(mod11, newx = X_test, s = 'lambda.min', gamma = 'gamma.min')


# 12. model
# on_off_comb_subject0_lasso_modlambda.min_gamma.min
f <- formula(paste('~ subject_id +', paste0(simple_int_acc_mag_features_names, collapse = ' + '), '+',
                   paste0(ts_features_names, collapse = ' + ')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

mod12 <- cv.glmnet(x = X, 
                   y = train_dat$on_off, 
                   weights = weights, 
                   family = 'gaussian',
                   penalty.factor = as.numeric(!grepl('subject', colnames(X)) | grepl(':', colnames(X))), 
                   alpha = 1, 
                   pmax = 200, 
                   relax = TRUE)

pred12 <- predict(mod12, newx = X_test, s = 'lambda.min', gamma = 'gamma.min')

#
# Combine models on subject level
#

cor(data.frame(pred1, pred2, pred3, pred4, pred5, pred6, pred7, pred8, pred9, pred10, pred11, pred12))
plot(data.frame(pred1, pred2, pred3, pred4, pred5, pred6, pred7, pred8, pred9, pred10, pred11, pred12), pch = '.')


# Mixed prediction model
best_model <- data.frame(subject_id = unique(train_dat$subject_id),
                         prediction_model = c(2, 7, 4, 6, 5, 7, 8, 7, 12, 9, 1, 11, 11, 10, 7))

mixed_pred_data <- data.frame(pred1, pred2, pred3, pred4, pred5, pred6, pred7, pred8, pred9, pred10, pred11, pred12,
                              subject_id = subset(test_predictions, subject_id %in% train_dat$subject_id)$subject_id)
mixed_pred_data$mixed <- as.numeric(mixed_pred_data[cbind(1:nrow(mixed_pred_data), best_model[match(mixed_pred_data$subject_id, best_model$subject_id), 'prediction_model'])])

# Censor predictions
mixed_pred_data$mixed[mixed_pred_data$mixed < 0] <- 0
mixed_pred_data$mixed[mixed_pred_data$mixed > 4] <- 4


cor(mixed_pred_data[, -13])
plot(mixed_pred_data[, c(1, 14)])

# Save results
test_predictions$on_off[test_predictions$subject_id %in% train_dat$subject_id] <- mixed_pred_data$mixed


#########################
# Dyskinesia best model #
#########################

# Training data
idx_train <- with(reg_data, !is.na(dyskinesia) & type != 'test')
train_dat <- subset(reg_data, idx_train)

train_dat <- subset(reg_data, !is.na(dyskinesia))

# Weights
tmp <- as.data.frame(table(train_dat$subject_id))
tmp[, 2] <- sqrt(tmp[, 2])
weights <- 1 / tmp[match(train_dat$subject_id, tmp[, 1]), 2]


# 1. model
mod1 <- lm(dyskinesia ~ subject_id, data = train_dat, weights = weights)
pred1 <- predict(mod1, newdata = subset(test_predictions, subject_id %in% train_dat$subject_id))

# 2. Model
# dyskinesia_simple_acc_mag_features_subject_xgb_mod
f <- formula(paste(' ~ subject_id +', paste0(simple_acc_mag_features_names, collapse = '+')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

# Fit model and predict on full data
mod2 <- xgboost(params = list(learning_rate = 0.3, nthread = 4), 
                X, weight = weights, label = train_dat$dyskinesia, nrounds = 10)

pred2 <- predict(mod2, newdata = X_test)

# 3. Model
# dyskinesia_simple_acc_mag_features_subject_rf_xgb_mod
f <- formula(paste(' ~ subject_id +', paste0(simple_acc_mag_features_names, collapse = '+')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

# Fit model and predict on full data
mod3 <- xgboost(params = list(colsample_bynode = 0.8, 
                              learning_rate = 1, 
                              max_depth = 5, 
                              num_parallel_tree = 500, 
                              subsample = 0.8,
                              nthread = 4),
                X, weight = weights, label = train_dat$dyskinesia, nrounds = 1)

pred3 <- predict(mod3, newdata = X_test)


# 4. Model
# dyskinesia_simple_int_acc_mag_features_subject_xgb_mod
f <- formula(paste(' ~ subject_id +', paste0(simple_int_acc_mag_features_names, collapse = '+')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

# Fit model and predict on full data
mod4 <- xgboost(params = list(learning_rate = 0.3, nthread = 4), 
                X, weight = weights, label = train_dat$dyskinesia, nrounds = 10)

pred4 <- predict(mod4, newdata = X_test)

# 5. Model
# dyskinesia_simple_int_acc_mag_features_subject_rf_xgb_mod
f <- formula(paste(' ~ subject_id +', paste0(simple_int_acc_mag_features_names, collapse = '+')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

# Fit model and predict on full data
mod5 <- xgboost(params = list(colsample_bynode = 0.8, 
                              learning_rate = 1, 
                              max_depth = 5, 
                              num_parallel_tree = 500, 
                              subsample = 0.8,
                              nthread = 4),
                X, weight = weights, label = train_dat$dyskinesia, nrounds = 1)

pred5 <- predict(mod5, newdata = X_test)

# 6. model
# dyskinesia_simple_acc_mag_features_subject0_ridge_mod
f <- formula(paste('~ subject_id +', paste0(simple_acc_mag_features_names, collapse = ' + ')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

mod6 <- cv.glmnet(x = X, 
                  y = train_dat$dyskinesia, 
                  weights = weights, 
                  family = 'gaussian',
                  penalty.factor = as.numeric(!grepl('subject', colnames(X)) | grepl(':', colnames(X))), 
                  alpha = 0)

pred6 <- predict(mod6, newx = X_test, s = 'lambda.min')


# Explore predictions
plot(data.frame(pred1, pred2, pred3, pred4, pred5, pred6), asp = 1)

cor(data.frame(pred1, pred2, pred3, pred4, pred5, pred6))

# Mixed prediction model
best_model <- data.frame(subject_id = unique(train_dat$subject_id),
                         prediction_model = c(3, 5, 1, 1, 6, 1, 3, 2, 1, 4, 1))

mixed_pred_data <- data.frame(pred1, pred2, pred3, pred4, pred5, pred6, subject_id = subset(test_predictions, subject_id %in% train_dat$subject_id)$subject_id)
mixed_pred_data$mixed <- as.numeric(mixed_pred_data[cbind(1:nrow(mixed_pred_data), best_model[match(mixed_pred_data$subject_id, best_model$subject_id), 'prediction_model'])])

# Censor predictions
mixed_pred_data$mixed[mixed_pred_data$mixed < 0] <- 0
mixed_pred_data$mixed[mixed_pred_data$mixed > 4] <- 4



plot(mixed_pred_data[, -7], asp = 1)

cor(mixed_pred_data[, -7])


# Save results
test_predictions$dyskinesia[test_predictions$subject_id %in% train_dat$subject_id] <- mixed_pred_data$mixed

######################
# Tremor best models #
######################

# Training data
idx_train <- with(reg_data, !is.na(tremor) & type != 'test')
train_dat <- subset(reg_data, idx_train)

train_dat <- subset(reg_data, !is.na(tremor))

# Weights
tmp <- as.data.frame(table(train_dat$subject_id))
tmp[, 2] <- sqrt(tmp[, 2])
weights <- 1 / tmp[match(train_dat$subject_id, tmp[, 1]), 2]


# 1. model
mod1 <- lm(tremor ~ subject_id, data = train_dat, weights = weights)
pred1 <- predict(mod1, newdata = subset(test_predictions, subject_id %in% train_dat$subject_id))

# 2. model
# tremor_all_per_subject_lasso_modlambda.1se_gamma.1se 

f <- formula(paste('~ -1 + ', paste0(unlist(sapply(paste0(used_feature_sets, '_names'), get)), collapse = ' + ')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

mod2 <- list()
pred2 <- rep(NA, nrow(X_test))

for (subj in unique(train_dat$subject_id)) {
  # Fit model and predict on subject data
  idx_subj <- with(train_dat, subject_id == subj)
  mod2[[subj]] <- cv.glmnet(x = X[idx_subj,], y = train_dat$tremor[idx_subj], 
                            family = 'gaussian', 
                            alpha = 1,
                            relax = TRUE)
  idx_subj <- with(subset(test_predictions, subject_id %in% train_dat$subject_id), subject_id == subj)
  
  pred2[idx_subj] <- predict(mod2[[subj]], newx = X_test[idx_subj, ], s = 'lambda.1se', gamma = 'gamma.1se')
}

# 3. Model
# tremor_comb_per_subject_lasso_modlambda.min_gamma.1se 
f <- formula(paste('~ -1 +', paste0(simple_int_acc_mag_features_names, collapse = ' + '), '+', 
                   paste0(ts_features_names, collapse = ' + ')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

mod3 <- list()
pred3 <- rep(NA, nrow(X_test))

for (subj in unique(train_dat$subject_id)) {
  # Fit model and predict on subject data
  idx_subj <- with(train_dat, subject_id == subj)
  mod3[[subj]] <- cv.glmnet(x = X[idx_subj,], y = train_dat$tremor[idx_subj], 
                            family = 'gaussian', 
                            alpha = 1, 
                            pmax = 200, 
                            relax = TRUE)
  idx_subj <- with(subset(test_predictions, subject_id %in% train_dat$subject_id), subject_id == subj)
  
  pred3[idx_subj] <- predict(mod3[[subj]], newx = X_test[idx_subj, ], s = 'lambda.min', gamma = 1)
}

# 4. model
# tremor_simple_int_acc_mag_features_subject0_pliable_0.5_mod 
f <- formula(paste('~ subject_id +', paste0(simple_int_acc_mag_features_names, collapse = ' + ')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

idx_keep <- apply(X, 2, sd) != 0
X <- X[, idx_keep]
X_test <- X_test[, idx_keep]

non_subj_cols <- !grepl('subject', colnames(X)) | grepl(':', colnames(X))

mx <- colMeans(X[, non_subj_cols])
sx <- sqrt(apply(X[, non_subj_cols], 2, var))
X[, non_subj_cols] <- scale(X[, non_subj_cols], mx, sx) 
X_test[, non_subj_cols] <- scale(X_test[, non_subj_cols], mx, sx) 


# Z design matrix
Z <- model.matrix(~ subject_id, train_dat)
Z_test <- model.matrix(~ subject_id, subset(test_predictions, subject_id %in% train_dat$subject_id))

idx_keep <- apply(Z, 2, sd) != 0
Z <- Z[, idx_keep]
Z_test <- Z_test[, idx_keep]

# Penalty factor
penalty.factor <- as.numeric(!grepl('subject', colnames(X)) | grepl(':', colnames(X)))

# Pliable lasso model
mod_pl <- pliable(x = X, 
                  z = Z, 
                  y = train_dat$tremor, 
                  alpha = 0.5,
                  w = weights,
                  family = 'gaussian', 
                  penalty.factor = penalty.factor, 
                  zlinear = FALSE,
                  maxinter = 2000)

mod_pl_cv <- cv.pliable(mod_pl, 
                        x = X, 
                        z = Z, 
                        y = train_dat$tremor, 
                        nfolds = 5,
                        verbose = FALSE)

pred4 <- predict(mod_pl, 
                 x = X_test, 
                 z = Z_test, 
                 lambda = mod_pl_cv$lambda.min)

# 5. model
# tremor_simple_acc_mag_features_subject0_ridge_mod 
f <- formula(paste('~ subject_id +', paste0(simple_acc_mag_features_names, collapse = ' + ')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

mod5 <- cv.glmnet(x = X, 
                  y = train_dat$tremor, 
                  weights = weights, 
                  family = 'gaussian',
                  penalty.factor = as.numeric(!grepl('subject', colnames(X)) | grepl(':', colnames(X))), 
                  alpha = 0)

pred5 <- predict(mod5, newx = X_test, s = 'lambda.min')


# 6. model
# tremor_comb_subject0_ridge_mod 
f <- formula(paste('~ subject_id +', paste0(simple_int_acc_mag_features_names, collapse = ' + '), '+',
                   paste0(ts_features_names, collapse = ' + ')))
# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

mod6 <- cv.glmnet(x = X, 
                  y = train_dat$tremor, 
                  weights = weights, 
                  family = 'gaussian',
                  penalty.factor = as.numeric(!grepl('subject', colnames(X)) | grepl(':', colnames(X))), 
                  alpha = 0)

pred6 <- predict(mod6, newx = X_test, s = 'lambda.min')

# 7. model
# tremor_PMA_features_subject0_ridge_mod 
f <- formula(paste('~ subject_id +', paste0(PMA_features_names, collapse = ' + ')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

mod7 <- cv.glmnet(x = X, 
                  y = train_dat$tremor, 
                  weights = weights, 
                  family = 'gaussian',
                  penalty.factor = as.numeric(!grepl('subject', colnames(X)) | grepl(':', colnames(X))), 
                  alpha = 0)

pred7 <- predict(mod7, newx = X_test, s = 'lambda.min')

# 8. model
# tremor_ts_int_features_subject0_lasso_modlambda.min_gamma.min 
f <- formula(paste('~ subject_id +', paste0(ts_int_features_names, collapse = ' + ')))

# Design matrix
X <- model.matrix(f, train_dat)
X_test <- model.matrix(f, subset(test_predictions, subject_id %in% train_dat$subject_id))

# Fit model and predict on full data
mod8 <- cv.glmnet(x = X, 
                  y = reg_data$tremor[idx_train], 
                  weights = weights, 
                  family = 'gaussian', 
                  penalty.factor = as.numeric(!grepl('subject', colnames(X))), 
                  alpha = 1,
                  relax = TRUE)


pred8 <- as.numeric(predict(mod8, newx = X_test, s = 'lambda.min', gamma = 'gamma.min'))

# Explore predictions
plot(data.frame(pred1, pred2, pred3, pred4, pred5, pred6, pred7, pred8), asp = 1)

cor(data.frame(pred1, pred2, pred3, pred4, pred5, pred6, pred7, pred8))

# Mixed prediction model
best_model <- data.frame(subject_id = unique(train_dat$subject_id),
                         prediction_model = c(1, 2, 1, 1, 3, 4, 5, 1, 1, 1, 8, 6, 7))

mixed_pred_data <- data.frame(pred1, pred2, pred3, pred4, pred5, pred6, pred7, pred8, subject_id = subset(test_predictions, subject_id %in% train_dat$subject_id)$subject_id)
mixed_pred_data$mixed <- as.numeric(mixed_pred_data[cbind(1:nrow(mixed_pred_data), best_model[match(mixed_pred_data$subject_id, best_model$subject_id), 'prediction_model'])])

# Censor predictions
mixed_pred_data$mixed[mixed_pred_data$mixed < 0] <- 0
mixed_pred_data$mixed[mixed_pred_data$mixed > 4] <- 4



plot(mixed_pred_data, asp = 1)

cor(mixed_pred_data[, -9])


# Save results
test_predictions$tremor[test_predictions$subject_id %in% train_dat$subject_id] <- mixed_pred_data$mixed

#
# Check covariance
#

cov(test_predictions[, c('on_off', 'dyskinesia', 'tremor')], use = 'pairwise.complete.obs')
cov(train_dat[, c('on_off', 'dyskinesia', 'tremor')], use = 'pairwise.complete.obs')

cor(test_predictions[, c('on_off', 'dyskinesia', 'tremor')], use = 'pairwise.complete.obs')
cor(train_dat[, c('on_off', 'dyskinesia', 'tremor')], use = 'pairwise.complete.obs')

#
# Test prediction
#

OnOff_Submission <- test_predictions[, c('measurement_id', 'on_off')]
names(OnOff_Submission)[2] <- 'prediction'

Dyskinesia_Submission <- test_predictions[, c('measurement_id', 'dyskinesia')]
names(Dyskinesia_Submission)[2] <- 'prediction'

Tremor_Submission <- test_predictions[, c('measurement_id', 'tremor')]
names(Tremor_Submission)[2] <- 'prediction'

