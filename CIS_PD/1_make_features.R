# Make data if it does not exist already
if (!exists('sensor_data')) source('0_load_data.R')

# Read packages
detach(package:plyr)
library(plyr)
library(tidyr)
library(mice)
library(flsa)
library(doParallel)
library(tsfeatures)
library(pracma)
library(PMA)


#
# Make time series features
#

features1d <- function(dat, variable) {
  dat %>% group_by(measurement_id) %>% mutate(XX = get(variable)) %>% summarize(
    # Time features
    mean_ts = mean(Timestamp), 
    max_ts = max(Timestamp),
    mean_diff_ts = mean(diff(Timestamp)),
    sd_ts = sd(Timestamp),
    
    obs = length(Timestamp),
    
    # 3D features
    mean = mean(XX),
    min = min(XX), 
    max = max(XX),
    abs_mean = mean(abs(XX)),
    sd = sd(XX),

    # Momenmts
    E1 = mean(XX),
    E2 = mean(XX^2),
    E3 = mean(XX^3),
    skew = mean(XX^3) / sqrt(mean(XX^2))^3,
    kuriosis = mean(XX^4) / mean(XX^2)^2,

    # Time-dependent moments
    Et1 = mean((Timestamp*XX)),
    Et2 = mean((Timestamp*XX)^2),
    Et3 = mean((Timestamp*XX)^3),
    skew_t = mean((Timestamp*XX)^3) / sqrt(mean((Timestamp*XX)^2))^3,
    kuriosis_t = mean((Timestamp*XX)^4) / mean((Timestamp*XX)^2)^2,

    # Differences
    meandiff = mean(diff(XX)), 
    meandiff_n = mean(diff(XX) / diff(Timestamp)), 
    meandiff2_n = mean(diff(XX, lag = 2) / diff(Timestamp, lag = 2)), 
    meandiff3_n = mean(diff(XX, lag = 3) / diff(Timestamp, lag = 3)), 
    meandiff5_n = mean(diff(XX, lag = 5) / diff(Timestamp, lag = 5)), 
    meandiff10_n = mean(diff(XX, lag = 10) / diff(Timestamp, lag = 10))) 
}

#
# Acceleration magnitude features
#

# Add acceleration magnitude to data
sensor_data$acc_mag <- with(sensor_data, X^2 + Y^2 + Z^2)

simple_acc_mag_features <- features1d(sensor_data, variable = 'acc_mag')
names(simple_acc_mag_features)[-1] <- paste0(names(simple_acc_mag_features)[-1], '_acc')

# Integration of acceleration magnitude over time

# Integrate acceleration magnitude with gravitation subtracted
sensor_data %<>% group_by(measurement_id) %>% mutate(int_acc_mag = cumtrapz(Timestamp, acc_mag - 1)) 

simple_int_acc_mag_features <- features1d(sensor_data, variable = 'int_acc_mag')
names(simple_int_acc_mag_features)[-1] <- paste0(names(simple_int_acc_mag_features)[-1], '_int_acc')

#
# Generated features
#

features <- c('acf_features', 
              'max_kl_shift', 
              # 'compengine', 
              'outlierinclude_mdrmd',
              'arch_stat', 
              'max_level_shift', 
              'ac_9', 
              # 'sampenc', # Very slow
              'crossing_points', 
              'max_var_shift', 
              # 'binarize_mean', 
              # 'sampen_first', # Very slow
              'entropy', 
              'nonlinearity', 
              'embed2_incircle', 
              'spreadrandomlocal_meantaul',
              'flat_spots', 
              'pacf_features', 
              'firstmin_ac', 
              'std1st_der',
              'heterogeneity', 
              'stability', 
              'firstzero_ac', 
              'trev_num',
              'holt_parameters', 
              'stl_features', 
              # 'fluctanal_prop_r1', # Very slow
              'walker_propcross',
              'hurst', 
              'unitroot_kpss', 
              'histogram_mode',
              # 'hw_parameters', 
              'unitroot_pp', 
              'localsimple_taures', 
              'lumpiness', 
              'motiftwo_entro3')

ts_features <- list()

for (id in unique(sensor_data$measurement_id)) {
  tmp <- with(subset(sensor_data, measurement_id == id), acc_mag)
  ts_features[[id]] <- as.matrix(tsfeatures(tmp, features = features))
}

ts_features <- do.call(rbind.data.frame, ts_features)
ts_features <- cbind(measurement_id = rownames(ts_features), ts_features)

# Integrated acceleration magnitude features
for (id in unique(sensor_data$measurement_id)) {
  tmp <- with(subset(sensor_data, measurement_id == id), int_acc_mag)
  # Special case for singular data
  if (id == 'dd261c90-aab6-40c3-9da7-1b8b320c7f7a') {
    tmp <- tmp + rnorm(length(tmp), sd = 0.001)
  }
  
  ts_int_features_l[[id]] <- as.matrix(tsfeatures(tmp, features = features))
}

ts_int_features <- do.call(rbind.data.frame, ts_int_features_l)
ts_int_features <- cbind(measurement_id = rownames(ts_int_features), ts_int_features)
names(ts_int_features)[-1] <- paste0(names(ts_int_features)[-1], '_int')

# UPDRS 3 - wide and imputed

UPDRS3_feat <- pivot_wider(UPDRS3, names_from = Visit, values_from = names(UPDRS3)[-(1:2)])
names(UPDRS3_feat) <- gsub(' ', '_', names(UPDRS3_feat))
names(UPDRS3_feat) <- gsub(':', '_', names(UPDRS3_feat))
UPDRS3_feat <- UPDRS3_feat[, colMeans(is.na(UPDRS3_feat)) < 1]

UPDRS3_imputed <- mice(UPDRS3_feat, seed = 1234, ridge = 0.001)
UPDRS3_feat <- complete(UPDRS3_imputed)
UPDRS3_feat <- UPDRS3_feat[, colSums(is.na(UPDRS3_feat)) == 0]

#
# Sparse aggregated features
#

SPC_feat_dat <- merge(merge(merge(simple_acc_mag_features, simple_int_acc_mag_features), ts_features), ts_int_features)
SPC_feat_dat <- SPC_feat_dat[, -(which(sapply(SPC_feat_dat[, -1], sd) == 0) + 1)]

# Scale
SPC_feat_dat[, -1] <- sapply(SPC_feat_dat[, -1], scale)

spc_cv <- SPC.cv(as.matrix(SPC_feat_dat[, -1]), sumabsvs = seq(1.1, floor(sqrt(ncol(SPC_feat_dat) - 1)), length = 20), niter = 20)

pmd_cv <- PMD.cv(as.matrix(SPC_feat_dat[, -1]), type = 'standard', niter = 20, sumabss = seq(0.09, 1, length = 50))


PMA_features <- data.frame(SPC_feat_dat[, 1],
                       SPC(as.matrix(SPC_feat_dat[, -1]), sumabsv = spc_cv$bestsumabsv, K = 20, orth = TRUE)$u,
                       PMD(as.matrix(SPC_feat_dat[, -1]), type = 'standard', sumabs = pmd_cv$bestsumabs, niter = 20, K = 20)$u)

names(PMA_features)[1] <- 'measurement_id'
names(PMA_features)[2:21] <- paste0('SPC', 1:20)
names(PMA_features)[22:41] <- paste0('PMD', 1:20)
