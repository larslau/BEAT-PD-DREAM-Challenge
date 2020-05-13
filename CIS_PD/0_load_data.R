library(readbulk)
library(data.table)

#
# Load training data (CIS-PD)
#

# Data directories
data_dir <- 'cis-pd/training_data/'
data_dir_test <- 'cis-pd/testing_data/'
data_dir_ancillary <- 'cis-pd/ancillary_data/'
data_dir_clinical <- 'cis-pd/clinical/'


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
names(sensor_data)[5] <- 'measurement_id'

# Build and save type data
type_dat <- sensor_data[, c('measurement_id', 'type')]
type_dat <- type_dat[!duplicated(type_dat$measurement_id), ]

# Clinical data
demographic <- read.csv(paste0(data_dir_clinical, 'CIS-PD_Demographics.csv'))
UPDRS124 <- read.csv(paste0(data_dir_clinical, 'CIS-PD_UPDRS_Part1_2_4.csv'))
UPDRS3 <- read.csv(paste0(data_dir_clinical, 'CIS-PD_UPDRS_Part3.csv'))

# Labels
labels <- read.csv(paste0(data_dir, 'CIS-PD_Training_Data_IDs_Labels.csv'))
labels <- rbind(labels, read.csv(paste0(data_dir, 'CIS-PD_Ancillary_Data_IDs_Labels.csv')))
labels <- merge(labels, read.csv(paste0(data_dir, 'cis-pd.CIS-PD_Test_Data_IDs.csv')), all = TRUE)
