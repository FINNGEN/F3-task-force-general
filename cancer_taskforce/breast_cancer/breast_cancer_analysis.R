# Analysis related to prelim check of breast cancer data for FG3 Cancer task force

library(bigrquery)
library(data.table)

# Source common functions
source("../../commons/operation_code_plots.R")

# BigQuery common params
projectid <- "finngen-refinery-dev"
db <- "sandbox_tools_r12"
options(gargle_oauth_cache = FALSE)
bq_auth(scopes = "https://www.googleapis.com/auth/bigquery")

minimum_extended <- tbl <- paste(projectid, db, "minimum_extended_r12_v1", sep = ".")
endpoint_longitudinal <- tbl <- paste(projectid, db, "endpoint_longitudinal_r12_v1", sep = ".")
detailed_longitudinal <- tbl <- paste(projectid, db, "finngen_r12_service_sector_detailed_longitudinal_v1", sep = ".")

# get individuals with breast cancer with their earliest diagnose age
sql <- paste(
    "SELECT FINNGENID, min(APPROX_EVENT_DAY) AS APPROX_EVENT_DAY, min(EVENT_AGE) as EVENT_AGE",
    "FROM", endpoint_longitudinal,
    "WHERE ENDPOINT LIKE 'C3_BREAST%'",
    "GROUP BY FINNGENID",
    "HAVING EVENT_AGE<80 AND APPROX_EVENT_DAY>='2005-01-01'"
)

breast_bqtable <- bq_project_query(projectid, sql)
breast_dt <- as.data.table(bq_table_download(breast_bqtable))

# get operations related to breast cancer
sql <- paste(
    "SELECT *",
    "FROM", detailed_longitudinal,
    "WHERE SOURCE like 'OPER_%' AND (CODE1 IN ('HA003','HA013','HA1AA','HA1AE','HA1CG','WF002','WF003','PJD52') OR CODE1 LIKE 'HAB%' OR CODE1 LIKE 'HAC%')"
)
oper_bqtable <- bq_project_query(projectid, sql)
oper_dt <- as.data.table(bq_table_download(oper_bqtable))

# Operations plot
setnames(x = breast_dt, old = c("APPROX_EVENT_DAY", "EVENT_AGE"), new = c("APPROX_DIAGNOSE_DAY", "DIAGNOSE_AGE"))
breast_oper <- merge(breast_dt, oper_dt, by = "FINNGENID")
breast_oper_filt <- breast_oper[EVENT_AGE > DIAGNOSE_AGE]

plot <- ggplot_geom_bar_operations(breast_oper_filt[, CODE1], lang = "fi")

#create view breast_min_event as select FINNGENID, min(APPROX_EVENT_DAY) AS APPROX_EVENT_DAY, min(EVENT_AGE) as EVENT_AGE from endpoint_long WHERE ENDPOINT LIKE 'C3_BREAST%' group by FINNGENID
#create view breast_min_event_filtered as select * from breast_min_event where EVENT_AGE<80 and APPROX_EVENT_DAY>='2005-01-01'


