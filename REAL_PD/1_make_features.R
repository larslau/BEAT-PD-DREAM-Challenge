# Make data if it does not exist already
if (!exists('sensor_data')) source('llra/REAL_PD/0_load_data.R')


# Read packages
detach(package:plyr)
detach(package:dplyr)
library(plyr)
library(dplyr)
library(tidyr)
library(mice)
library(flsa)
library(doParallel)
library(tsfeatures)
library(pracma)
library(PMA)


###########################################
# Load stored data instead of recomputing #
###########################################

used_feature_sets <- c('simple_acc_mag_features',
                       'simple_int_acc_mag_features',
                       'simple_gyro_mag_features',
                       'simple_int_gyro_mag_features',
                       'ts_features', 
                       'ts_int_features',
                       'ts_gyro_features',
                       'ts_gyro_int_features',
                       'PMA_features')


# Make data if it does not exist already
UPDRS3_feat <- read.csv('REAL_PD/features/UPDRS3_features.csv')

# Load features
for (feat in used_feature_sets) {
  assign(feat, read.csv(paste0('REAL_PD/features/', feat, '.csv')))
}

# 
# #
# # Make time series features
# #
# 
# features1d <- function(dat, variable) {
#   dat %>% group_by(measurement_id) %>% mutate(XX = get(variable)) %>% summarize(
#     # Time features
#     mean_ts = mean(Timestamp), 
#     max_ts = max(Timestamp),
#     mean_diff_ts = mean(diff(Timestamp)),
#     sd_ts = sd(Timestamp),
#     
#     obs = length(Timestamp),
#     
#     # 3D features
#     mean = mean(XX),
#     min = min(XX), 
#     max = max(XX),
#     abs_mean = mean(abs(XX)),
#     sd = sd(XX),
# 
#     # Momenmts
#     E1 = mean(XX),
#     E2 = mean(XX^2),
#     E3 = mean(XX^3),
#     skew = mean(XX^3) / sqrt(mean(XX^2))^3,
#     kuriosis = mean(XX^4) / mean(XX^2)^2,
# 
#     # Time-dependent moments
#     Et1 = mean((Timestamp*XX)),
#     Et2 = mean((Timestamp*XX)^2),
#     Et3 = mean((Timestamp*XX)^3),
#     skew_t = mean((Timestamp*XX)^3) / sqrt(mean((Timestamp*XX)^2))^3,
#     kuriosis_t = mean((Timestamp*XX)^4) / mean((Timestamp*XX)^2)^2,
# 
#     # Differences
#     meandiff = mean(diff(XX)), 
#     meandiff_n = mean(diff(XX) / ifelse(diff(Timestamp) == 0, 0.005, diff(Timestamp))), 
#     meandiff2_n = mean(diff(XX, lag = 2) / ifelse(diff(Timestamp, lag = 2) == 0, 0.005, diff(Timestamp, lag = 2))), 
#     meandiff3_n = mean(diff(XX, lag = 3) / ifelse(diff(Timestamp, lag = 3) == 0, 0.005, diff(Timestamp, lag = 3))), 
#     meandiff5_n = mean(diff(XX, lag = 5) / ifelse(diff(Timestamp, lag = 5) == 0, 0.005, diff(Timestamp, lag = 5))), 
#     meandiff10_n = mean(diff(XX, lag = 10) / ifelse(diff(Timestamp, lag = 10) == 0, 0.005, diff(Timestamp, lag = 10)))) 
# }
# 
# #
# # Acceleration magnitude features
# #
# 
# # Add acceleration magnitude to data
# sensor_data$acc_mag <- with(sensor_data, X^2 + Y^2 + Z^2)
# 
# simple_acc_mag_features <- features1d(sensor_data, variable = 'acc_mag')
# names(simple_acc_mag_features)[-1] <- paste0(names(simple_acc_mag_features)[-1], '_acc')
# 
# # Integration of acceleration magnitude over time
# 
# # Integrate acceleration magnitude with gravitation subtracted
# sensor_data %<>% group_by(measurement_id) %>% mutate(int_acc_mag = cumtrapz(Timestamp, acc_mag - 1)) 
# # normalize to avoid numerical problems
# sensor_data$int_acc_mag <- sensor_data$int_acc_mag / diff(range(sensor_data$int_acc_mag))
# 
# simple_int_acc_mag_features <- features1d(sensor_data, variable = 'int_acc_mag')
# names(simple_int_acc_mag_features)[-1] <- paste0(names(simple_int_acc_mag_features)[-1], '_int_acc')
# 
# #
# # Gyroscope magnitude features 
# #
# 
# # Add acceleration magnitude to data
# gyro_data$gyro_mag <- with(gyro_data, sqrt(X^2 + Y^2 + Z^2))
# 
# simple_gyro_mag_features <- features1d(gyro_data, variable = 'gyro_mag')
# names(simple_gyro_mag_features)[-1] <- paste0(names(simple_gyro_mag_features)[-1], '_gyro')
# 
# # Integration of acceleration magnitude over time
# 
# # Integrate acceleration magnitude with gravitation subtracted
# gyro_data %<>% group_by(measurement_id) %>% mutate(int_gyro_mag = cumtrapz(Timestamp, gyro_mag - 1)) 
# # normalize to avoid numerical problems
# gyro_data$int_gyro_mag <- gyro_data$int_gyro_mag / diff(range(gyro_data$int_gyro_mag))
# 
# simple_int_gyro_mag_features <- features1d(gyro_data, variable = 'int_gyro_mag')
# names(simple_int_gyro_mag_features)[-1] <- paste0(names(simple_int_gyro_mag_features)[-1], '_int_gyro')
# 
# #
# # Acceleration magnitude features - phone 
# #
# 
# # Add acceleration magnitude to data
# phone_data$acc_mag <- with(phone_data, X^2 + Y^2 + Z^2)
# names(phone_data)[6] <- 'measurement_id'
# 
# simple_acc_mag_p_features <- features1d(phone_data, variable = 'acc_mag')
# names(simple_acc_mag_p_features)[1] <- 'measurement_id_original'
# names(simple_acc_mag_p_features)[-1] <- paste0(names(simple_acc_mag_p_features)[-1], '_acc_p')
# 
# # Integration of acceleration magnitude over time
# 
# # Integrate acceleration magnitude with gravitation subtracted
# phone_data %<>% group_by(measurement_id) %>% mutate(int_acc_mag = cumtrapz(Timestamp, acc_mag - 1)) 
# # normalize to avoid numerical problems
# phone_data$int_acc_mag <- phone_data$int_acc_mag / diff(range(phone_data$int_acc_mag))
# 
# simple_int_acc_mag_p_features <- features1d(phone_data, variable = 'int_acc_mag')
# names(simple_int_acc_mag_p_features)[1] <- 'measurement_id_original'
# names(simple_int_acc_mag_p_features)[-1] <- paste0(names(simple_int_acc_mag_p_features)[-1], '_int_acc_p')
# 
# #
# # Generated features
# #
# 
# features <- c('acf_features', 
#               'max_kl_shift', 
#               # 'compengine', 
#               'outlierinclude_mdrmd',
#               'arch_stat', 
#               'max_level_shift', 
#               'ac_9', 
#               # 'sampenc', # Very slow
#               'crossing_points', 
#               'max_var_shift', 
#               # 'binarize_mean', 
#               # 'sampen_first', # Very slow
#               'entropy', 
#               'nonlinearity', 
#               'embed2_incircle', 
#               'spreadrandomlocal_meantaul',
#               'flat_spots', 
#               'pacf_features', 
#               'firstmin_ac', 
#               'std1st_der',
#               'heterogeneity', 
#               'stability', 
#               'firstzero_ac', 
#               'trev_num',
#               'holt_parameters', 
#               'stl_features', 
#               # 'fluctanal_prop_r1', # Very slow
#               'walker_propcross',
#               'hurst', 
#               'unitroot_kpss', 
#               'histogram_mode',
#               # 'hw_parameters', 
#               'unitroot_pp', 
#               'localsimple_taures', 
#               'lumpiness', 
#               'motiftwo_entro3')
# 
# ts_features <- list()
# 
# for (id in unique(sensor_data$measurement_id)) {
#   tmp <- with(subset(sensor_data, measurement_id == id), acc_mag)
#   ts_features[[id]] <- as.matrix(tsfeatures(tmp, features = features))
# }
# 
# ts_features <- do.call(rbind.data.frame, ts_features)
# ts_features <- cbind(measurement_id = rownames(ts_features), ts_features)
# 
# # Integrated acceleration magnitude features
# for (id in unique(sensor_data$measurement_id)) {
#   tmp <- with(subset(sensor_data, measurement_id == id), int_acc_mag)
#   ts_int_features_l[[id]] <- as.matrix(tsfeatures(tmp, features = features))
# }
# 
# ts_int_features <- do.call(rbind.data.frame, ts_int_features_l)
# ts_int_features <- cbind(measurement_id = rownames(ts_int_features), ts_int_features)
# names(ts_int_features)[-1] <- paste0(names(ts_int_features)[-1], '_int')
# 
# #
# # Gyroscope features
# #
# 
# ts_gyro_features <- list()
# 
# for (id in unique(gyro_data$measurement_id)) {
#   tmp <- with(subset(gyro_data, measurement_id == id), gyro_mag)
#   ts_gyro_features[[id]] <- as.matrix(tsfeatures(tmp, features = features))
# }
# 
# ts_gyro_features <- do.call(rbind.data.frame, ts_gyro_features)
# ts_gyro_features <- cbind(measurement_id = rownames(ts_gyro_features), ts_gyro_features)
# names(ts_gyro_features)[-1] <- paste0(names(ts_gyro_features)[-1], '_gyro')
# 
# # Integrated gyroscope magnitude features
# 
# for (id in unique(gyro_data$measurement_id)) {
#   tmp <- with(subset(gyro_data, measurement_id == id), int_gyro_mag)
#   ts_gyro_int_features_l[[id]] <- as.matrix(tsfeatures(tmp, features = features))
# }
# 
# ts_gyro_int_features <- do.call(rbind.data.frame, ts_gyro_int_features_l)
# ts_gyro_int_features <- cbind(measurement_id = rownames(ts_gyro_int_features), ts_gyro_int_features)
# names(ts_gyro_int_features)[-1] <- paste0(names(ts_gyro_int_features)[-1], '_int_gyro')
# 
# #
# # Phone accelerometer features
# #
# 
# ts_phone_features <- list()
# 
# for (id in unique(phone_data$measurement_id)) {
#   tmp <- with(subset(phone_data, measurement_id == id), acc_mag)
#   ts_phone_features[[id]] <- as.matrix(tsfeatures(tmp, features = features))
# }
# 
# ts_phone_features <- do.call(rbind.data.frame, ts_phone_features)
# ts_phone_features <- cbind(measurement_id = rownames(ts_phone_features), ts_phone_features)
# names(ts_phone_features)[1] <- 'measurement_id_original'
# names(ts_phone_features)[-1] <- paste0(names(ts_phone_features)[-1], '_phone')
# 
# 
# # Integrated gyroscope magnitude features
# for (id in unique(phone_data$measurement_id)) {
#   tmp <- with(subset(phone_data, measurement_id == id), int_gyro_mag)
#   ts_phone_int_features_l[[id]] <- as.matrix(tsfeatures(tmp, features = features))
# }
# 
# ts_phone_int_features <- do.call(rbind.data.frame, ts_phone_int_features_l)
# ts_phone_int_features <- cbind(measurement_id = rownames(ts_phone_int_features), ts_phone_int_features)
# names(ts_phone_int_features)[1] <- 'measurement_id_original'
# names(ts_phone_int_features)[-1] <- paste0(names(ts_phone_int_features)[-1], '_phone_int')
# 
# # UPDRS 3 - wide
# 
# UPDRS3_feat <- pivot_wider(UPDRS3, names_from = ParticipantState, values_from = names(UPDRS3)[-(1:2)])
# names(UPDRS3_feat) <- gsub(' ', '_', names(UPDRS3_feat))
# names(UPDRS3_feat) <- gsub(':', '_', names(UPDRS3_feat))
# UPDRS3_feat <- UPDRS3_feat[, colMeans(is.na(UPDRS3_feat)) <= 0]
# 
# 
# #
# # Sparse aggregated features
# #
# 
# used_feature_sets <- c('simple_acc_mag_features',
#                        'simple_int_acc_mag_features',
#                        'simple_gyro_mag_features',
#                        'simple_int_gyro_mag_features',
#                        'simple_acc_mag_p_features',
#                        'simple_int_acc_mag_p_features',
#                        'ts_features',
#                        'ts_int_features',
#                        'ts_gyro_features',
#                        'ts_gyro_int_features',
#                        'ts_phone_features',
#                        'ts_phone_int_features')
# 
# SPC_feat_dat <- merge(merge(merge(labels, demographic, all = TRUE), UPDRS3_feat, all = TRUE), UPDRS124, all = TRUE)
# for (feat in used_feature_sets) {
#   SPC_feat_dat <- merge(SPC_feat_dat, get(feat), all = TRUE)
# }
# 
# # Remove factors and outcomes
# SPC_feat_dat <- SPC_feat_dat[, -which(names(SPC_feat_dat) %in% c('subject_id', 'Gender', 'on_off', 'dyskinesia', 'tremor'))]
# 
# SPC_feat_dat <- SPC_feat_dat[, -(which(sapply(SPC_feat_dat[, -(1:2)], sd, na.rm = TRUE) == 0) + 2)]
# 
# # Scale
# SPC_feat_dat[, -(1:2)] <- sapply(SPC_feat_dat[, -(1:2)], scale)
# 
# spc_cv <- SPC.cv(as.matrix(SPC_feat_dat[, -(1:2)]), sumabsvs = seq(1.1, floor(sqrt(ncol(SPC_feat_dat) - 1)), length = 20), niter = 20)
# 
# pmd_cv <- PMD.cv(as.matrix(SPC_feat_dat[, -(1:2)]), type = 'standard', niter = 20, sumabss = seq(0.09, 1, length = 50))
# 
# 
# PMA_features <- data.frame(SPC_feat_dat[, (1:2)],
#                        SPC(as.matrix(SPC_feat_dat[, -(1:2)]), sumabsv = spc_cv$bestsumabsv, K = 20, orth = TRUE)$u,
#                        PMD(as.matrix(SPC_feat_dat[, -(1:2)]), type = 'standard', sumabs = pmd_cv$bestsumabs, niter = 20, K = 20)$u)
# 
# names(PMA_features)[1] <- 'measurement_id_original'
# names(PMA_features)[2] <- 'measurement_id'
# names(PMA_features)[3:22] <- paste0('SPC', 1:20)
# names(PMA_features)[23:42] <- paste0('PMD', 1:20)
