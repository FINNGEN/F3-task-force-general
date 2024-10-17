#' @title Plasma Sample Data Basic Query
#' @description This script combines plasma availability information with endpoint data and performs basic queries.
#' @details The script performs the following tasks:
#' - Authenticates and sets up the BigQuery project.
#' - Defines a function to calculate approximate event dates given birth dates and ages at events.
#' - Gathers sample annotations from a BigQuery table.
#' - Reads plasma sample availability data from multiple files and combines them.
#' - Processes and cleans the combined plasma sample data.
#' - Retrieves phenotype data for specific endpoints from a BigQuery table.
#' - Adds additional kidney disease sample information.
#' - Merges phenotype data with plasma sample data.
#' - Performs an example query to get unique individuals who meet specific criteria related to plasma samples and endpoints.
#' 
#' @section General Settings:
#' - `projectid`: The BigQuery project ID (Sandbox name; Check the URL in your web browser).
#' - `options(gargle_oauth_cache = FALSE)`: Disables OAuth cache.
#' - `bq_auth(scopes = "https://www.googleapis.com/auth/bigquery")`: Authenticates BigQuery access.
#' 
#' @section Functions:
#' - `calc_approx_date(date, event_age)`: Calculates approximate event dates given birth dates and ages at events.
#' 
#' @section Data Gathering:
#' - Retrieves sample annotations from the `finngen-production-library.sandbox_tools_r12.minimum_extended_r12_v1` table.
#' - Reads plasma sample availability data from multiple files and combines them.
#' - Processes and cleans the combined plasma sample data.
#' 
#' @section Phenotype Data:
#' - Retrieves phenotype data for specific endpoints from the `finngen-production-library.sandbox_tools_r12.endpoint_cohorts_r12_v1` table.
#' - Adds additional kidney disease sample information from a local file.
#' 
#' @section Data Merging:
#' - Merges phenotype data with plasma sample data.
#' 
#' @section Example Query:
#' - Retrieves unique individuals who:
#'   - Have the endpoint `ILD_ENDPOINTS`.
#'   - Have plasma samples available.
#'   - Have plasma samples frozen within 8 hours of collection.
#'   - Have plasma samples taken before the first diagnosis.
#'   - Have plasma samples taken after 2013.
#' 
#' @import bigrquery
#' @import data.table
#' @import lubridate
#' @import ggplot2

library(bigrquery)
library(data.table)
library(lubridate)
library(ggplot2)

# General settings

# projectid is your sandbox name (check the URL in your web browser)
projectid <- "fg-production-sandbox-6"
options(gargle_oauth_cache = FALSE)
bq_auth(scopes = "https://www.googleapis.com/auth/bigquery")

#' Calculate Approximate Event Date
#'
#' This function calculates the approximate event date given a birth date and age at the event.
#'
#' @param date A date object or a character string representing the birth date.
#' @param event_age A numeric value representing the age at the event.
#'
#' @return A date object representing the approximate event date.
#'
#' @examples
#' calc_approx_date("2000-01-01", 20.5)
#' calc_approx_date(as.Date("2000-01-01"), 20.5)
#'
#' @importFrom lubridate ymd years days
#' @export
calc_approx_date <- function(date, event_age) {
  if (is.character(date)) {
    date <- ymd(date)
  }
  int_p <- floor(event_age)
  fract_p <- event_age - int_p
  approx_days <- round(fract_p * 365.25)
  approx_event_date <- date + years(int_p) + days(approx_days)
  return(approx_event_date)
}

# Gather some sample annotations

sql <- paste0(
  "SELECT * ",
  "FROM finngen-production-library.",
  "sandbox_tools_r12.",
  "minimum_extended_r12_v1 ",
  "AS A ",
  "LEFT JOIN (",
  "SELECT IID, batch ",
  "FROM finngen-production-library.",
  "sandbox_tools_r12.",
  "covariates_r12_v1)",
  "AS B ",
  "ON A.FINNGENID = B.IID"
)
anno <- as.data.table(bq_table_download(bq_project_query(projectid, sql)))
anno[, IN_DF12_ANALYSIS := ! is.na(batch)]

# Read in plasma sample availability data

F10_1 <- fread("/finngen/library-red/task_force_data/plasma_availability/F10/data/F10_1.txt")
F10_2 <- fread("/finngen/library-red/task_force_data/plasma_availability/F10/data/F10_2.txt")
F06 <- fread("/finngen/library-red/task_force_data/plasma_availability/F06_F07_v2/data/F06_v2.txt")
F07 <- fread("/finngen/library-red/task_force_data/plasma_availability/F06_F07_v2/data/F07_v2.txt")
BS <- fread("/finngen/library-red/task_force_data/plasma_availability/blood_service_plasma_sample_availability/data/blood_service_plasma_sample_availability.txt")

setnames(BS, "BIOBANK", "BIOBANK_PLASMA")

RR <- rbind(F06, F07, F10_1, F10_2, BS, fill = T)
RR[, COHORT_FINNGENID := NULL] # COHORT_FINNGENID not set for all data so remove and add later from annotation
RR <- RR[!is.na(FINNGENID)]
RR <- unique(RR)
RR[, FINNGENID_N_PLASMA_SAMPLES := .N, by = .(FINNGENID)]
RR[, SERIAL_SAMPLES := FINNGENID_N_PLASMA_SAMPLES > 1]
RR[, HOURS_FROM_COLLECTION_TO_FREEZING := round(time_length(difftime(APPROX_TIMESTAMP_FREEZING, APPROX_TIMESTAMP_COLLECTION), "hours"), 2)]
RR[, HOURS_FROM_COLLECTION_TO_FREEZING_BIN := cut(HOURS_FROM_COLLECTION_TO_FREEZING, breaks = c(0, 4, 6, 8, Inf), labels = c("<4h", "<6h", "<8h", ">8h"))]
RR <- merge(RR, anno, by = "FINNGENID", all = T)
RR[is.na(AVAILABILITY), ":="(AVAILABILITY = "NO", FINNGENID_N_PLASMA_SAMPLES = 0)]
RR[, BIOBANK := gsub("^THL BIOBANK.*", "THL BIOBANK", COHORT)]
RR[, BIOBANK := gsub("^ARCTIC BIOBANK.*", "ARCTIC BIOBANK", BIOBANK)]

# Get phenotype data

endpoints <- strsplit(c(
"K11_KELAIBD
C3_CANCER_WIDE
G6_MS
G6_AD_WIDE
G6_PARKINSON
G6_ALS
K11_FIBROCHIRLIV
H7_AMD
ILD_ENDPOINTS
I9_HEARTFAIL
M13_ANKYLOSPON
M13_SLE
M13_RHEUMA
M13_PSORIARTH
M13_SACROILIITIS
M13_SYSTSLCE
N14_RENFAIL
J10_ASTHMA
J10_COPD
N14_ACUTERENFAIL
N14_CHRONKIDNEYDIS
DM_RETINOPATHY
L12_ATOPIC"), split = "\n")[[1]]

sql <- paste0(
  "SELECT FINNGENID, ENDPOINT, AGE, APPROX_EVENT_DAY ",
  "FROM finngen-production-library.",
  "sandbox_tools_r12.",
  "endpoint_cohorts_r12_v1 ",
  "WHERE ENDPOINT IN (",
  gsub(",$", ")", paste(paste0("'", endpoints, "',"), collapse = " ")),
  " AND CONTROL_CASE_EXCL=1"
)

DFX <- as.data.table(bq_table_download(bq_project_query(projectid, sql)))
DFX[, APPROX_EVENT_DAY := as_datetime(APPROX_EVENT_DAY, tz = "Europe/Helsinki")]
setnames(DFX, old = "AGE", new = "EVENT_AGE")

# Add Kidney disease TF sample info
egfrdecline <- fread("/finngen/red/rwalters/ckd/EGFRDECLINE25.gate.pheno.txt")
egfrdecline[, ENDPOINT := "CKD_EGFRDECLINE25"]
setnames(x = egfrdecline, old = "FID", new = "FINNGENID")
egfrdecline[, EVENT_AGE := as.numeric(NA)]
egfrdecline[, APPROX_EVENT_DAY := as.POSIXct(NA)]

DFX <- rbind(DFX, egfrdecline[, .(FINNGENID, ENDPOINT, EVENT_AGE, APPROX_EVENT_DAY)])

# Merge phenotype data with plasma sample data

DATA <- merge(x = DFX,
              y = RR,
              by = "FINNGENID",
              all.x = T)

# Example query: Get all unique individuals who
# 1) have endpoint ILD_ENDPOINTS (ENDPOINT=="ILD_ENDPOINTS"), and
# 2) have plasma sample (AVAILABILITY=="YES"), and
# 3) plasma sample is frozen within 8h of collection (HOURS_FROM_COLLECTION_TO_FREEZING < 8), and
# 4) plasma sample is taken before first diagnosis (APPROX_TIMESTAMP_COLLECTION < APPROX_EVENT_DAY), and
# 5) plasma sample is taken after 2013

data_tmp <- DATA[
  ENDPOINT == "ILD_ENDPOINTS"
  & AVAILABILITY == "YES"
  & HOURS_FROM_COLLECTION_TO_FREEZING < 8
  & APPROX_TIMESTAMP_COLLECTION < APPROX_EVENT_DAY
  & year(APPROX_TIMESTAMP_COLLECTION) > 2013
]
uniq_fgids <- unique(data_tmp[, .(FINNGENID)])
