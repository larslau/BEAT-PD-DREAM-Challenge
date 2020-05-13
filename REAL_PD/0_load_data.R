library(readbulk)
library(data.table)

#
# Load training data (REAL-PD)
#

data_dir <- 'real-pd/training_data_updated/smartwatch_accelerometer/'
data_dir_test <- 'real-pd/testing_data/smartwatch_accelerometer/'
data_dir_ancillary <- 'real-pd/ancillary_data_updated/smartwatch_accelerometer/'
data_dir_clinical <- 'real-pd/clinical/'


# Sensor data
training_data <- read_bulk(directory = data_dir)
ancillary_data <- read_bulk(directory = data_dir_ancillary)
test_data <- read_bulk(directory = data_dir_test)

# Merge training, ancillary and test
training_data$type <- 'train'
ancillary_data$type <- 'ancillary'
test_data$type <- 'test'

sensor_data <- rbind(training_data, ancillary_data, test_data)

# Remove .csv from id
sensor_data$File <- gsub('.csv', '', sensor_data$File)
names(sensor_data)[6] <- 'measurement_id'

names(sensor_data)[1] <- 'Timestamp'
names(sensor_data)[3:5] <- c('X', 'Y', 'Z')

# New measurement id
sensor_data$measurement_id_original <- sensor_data$measurement_id
sensor_data$measurement_id <- with(sensor_data, paste0(measurement_id, '_', device_id))

# Remove samples with less than 5000 measurements
sensor_data <- subset(sensor_data, measurement_id %in% names(which(table(measurement_id) >= 5000)))

# Build and save type data
type_dat <- sensor_data[, c('measurement_id', 'measurement_id_original', 'device_id', 'type')]
type_dat <- type_dat[!duplicated(type_dat$measurement_id), ]

# Clinical data
demographic <- read.csv(paste0(data_dir_clinical, 'REAL-PD_Demographics.csv'))
UPDRS124 <- read.csv(paste0(data_dir_clinical, 'clinical_data/REAL-PD_UPDRS_Part1_2_4.csv'))
UPDRS3 <- read.csv(paste0(data_dir_clinical, 'clinical_data/REAL-PD_UPDRS_Part3.csv'))

# Smartphone metadata
smartphone_metadata <- read.csv(paste0(data_dir, 'clinical_data/REAL-PD_Smartphone_Metadata.csv'))


# Labels
labels <- read.csv(paste0(data_dir, 'REAL-PD_Training_Data_IDs_Labels.csv'))
labels <- rbind(labels, read.csv(paste0(data_dir, 'REAL-PD_Ancillary_Data_IDs_Labels.csv')))
labels <- merge(labels, read.csv(paste0(data_dir, 'real-pd.REAL-PD_Test_Data_IDs.csv')), all = TRUE)
names(labels)[1] <- 'measurement_id_original'
labels <- merge(labels, type_dat[, c('measurement_id', 'measurement_id_original')], all = TRUE)

#
# Gyroscope data
#

data_dir_gyro <- 'real-pd/training_data_updated/smartwatch_gyroscope/'
data_dir_test_gyro <- 'real-pd/testing_data/smartwatch_gyroscope/'
data_dir_ancillary_gyro <- 'real-pd/ancillary_data_updated/smartwatch_gyroscope/'


# Sensor data
training_data <- read_bulk(directory = data_dir_gyro)
ancillary_data <- read_bulk(directory = data_dir_ancillary_gyro)
test_data <- read_bulk(directory = data_dir_test_gyro)

# Merge training, ancillary and test
training_data$type <- 'train'
ancillary_data$type <- 'ancillary'
test_data$type <- 'test'

gyro_data <- rbind(training_data, ancillary_data, test_data)

# Remove .csv from id
gyro_data$File <- gsub('.csv', '', gyro_data$File)
names(gyro_data)[6] <- 'measurement_id'

names(gyro_data)[1] <- 'Timestamp'
names(gyro_data)[3:5] <- c('X', 'Y', 'Z')

# New measurement id
gyro_data$measurement_id_original <- gyro_data$measurement_id
gyro_data$measurement_id <- with(gyro_data, paste0(measurement_id, '_', device_id))

# Remove samples with less than 5000 measurements
gyro_data <- subset(gyro_data, measurement_id %in% names(which(table(measurement_id) >= 5000)))

# Build and save type data
type_dat_g <- gyro_data[, c('measurement_id', 'measurement_id_original', 'device_id', 'type')]
type_dat_g <- type_dat_g[!duplicated(type_dat_g$measurement_id), ]

# Check
identical(dim(type_dat_g), dim(type_dat))
table(type_dat_g$measurement_id %in% type_dat$measurement_id)

#
# Smartphone data
#


data_dir_phone <- 'real-pd/training_data/smartphone_accelerometer/'
data_dir_test_phone <- 'real-pd/testing_data/smartphone_accelerometer/'
data_dir_ancillary_phone <- 'real-pd/ancillary_data/smartphone_accelerometer/'


# Sensor data
training_data <- read_bulk(directory = data_dir_phone)
ancillary_data <- read_bulk(directory = data_dir_ancillary_phone)
test_data <- read_bulk(directory = data_dir_test_phone)

# Merge training, ancillary and test
training_data$type <- 'train'
ancillary_data$type <- 'ancillary'
test_data$type <- 'test'

phone_data <- rbind(training_data, ancillary_data, test_data)

# Remove .csv from id
phone_data$File <- gsub('.csv', '', phone_data$File)
names(phone_data)[5] <- 'measurement_id'

names(phone_data)[1] <- 'Timestamp'
names(phone_data)[2:4] <- c('X', 'Y', 'Z')

# New measurement id
phone_data$measurement_id_original <- phone_data$measurement_id
phone_data$measurement_id <- NULL

# Remove samples with less than 5000 measurements
phone_data <- subset(phone_data, measurement_id_original %in% names(which(table(measurement_id_original) >= 5000)))


# Build and save type data
type_dat_p <- phone_data[, c('measurement_id_original', 'type')]
type_dat_p <- type_dat_p[!duplicated(type_dat_p$measurement_id), ]

# Check
dim(type_dat_p)
dim(type_dat)

table(type_dat_p$measurement_id_original %in% type_dat$measurement_id_original)

type_dat <- merge(type_dat, type_dat_p, all = TRUE)
