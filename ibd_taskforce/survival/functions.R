library(survival)
library(ggplot2)
library(survminer)
library(tidyverse)
library(lubridate)

#' Create survival analysis data from FINNGEN endpoint status and approximate event days
#' This function expects that the samples have a condition that might lead to an event, and a possible censoring event
#' Conditions are filtered by time to be between start and end dates
#'
#' @param data A dataframe 
#' @param condition Condition column that starts the followup period
#' @param event Event column
#' @param start_date Minimum date for condition
#' @param end_date Maximum date for condition
#' @param followup_period Followup period for events in years, by default 5. E.g. value of 1 means a 1-year followup from condition.
#' @param censor Censoring column, by default DEATH
#' @param additional_columns A vector of columns to extract into the survival data, to be included as covariates
#' @param stratify_into_n Stratify the time to condition from start date to n equal values. Useful for kaplan-meyer plots 
#' @returns A data frame with columns 'status', 'time', 'date_modifier', 'date_interval', and possible additional columns.
survival_analysis_helper <- function(data,
                                    condition,
                                    event,
                                    start_date,
                                    end_date,
                                    followup_period=5,
                                    censor="DEATH",
                                    additional_columns = c(),
                                    stratify_into_n = 4)
    {
    ## helper function to filter and transform diagnoses and diagnosis dates into survival time and status
    ## columns for condition and event dates are taken as column name + "_APPROX_EVENT_DAY"
    condition_date = paste0(condition,"_APPROX_EVENT_DAY",sep="")
    event_date = paste0(event,"_APPROX_EVENT_DAY",sep="")
    censor_date = paste0(censor,"_APPROX_EVENT_DAY",sep="")
    # filter data to correct samples
    cohort <- data %>% filter(get({{condition}}) == 1 &
                              get({{condition_date}}) >= as.Date(start_date) &
                              get({{condition_date}}) <= as.Date(end_date) &
                              (is.na(get({{event_date}})) | get({{event_date}}) >= get({{condition_date}}) ) &
                              (is.na(get({{censor_date}}))|get({{censor_date}}) > get({{condition_date}})))
    cohort <- cohort %>% rowwise() %>% mutate(status = if_else((!is.na(get({{event_date}})))&
                                                               (get({{event_date}}) < get({{condition_date}}) %m+%  years(followup_period))&
                                                               (is.na(get({{censor_date}}))|get({{event_date}}) < get({{censor_date}}) ),1,0),
                                              end_date = min(get({{censor_date}}),get({{condition_date}}) %m+% years(followup_period),get({{event_date}}) ,na.rm=T),
                                              time = difftime(end_date,get({{condition_date}}),units="days"),
                                              condition_date = get({{condition_date}})) %>% ungroup()
    first_condition = min(cohort[[condition_date]], na.rm = T)
    # create breaks to stratify the date modifier into {{stratify_into_n}} even buckets. 
    survival_data <-cohort %>% mutate(date_modifier = as.numeric( difftime(get({{condition_date}}),first_condition, units="days") )/365.25 ) %>%
                    select(status,time,date_modifier,additional_columns) %>% 
                    mutate(date_interval = cut(date_modifier,breaks=stratify_into_n))
    return(survival_data)
}

#' Collapse multiple columns to condition column, taking the first as the condition date
#'
#' @param data Dataframe
#' @param columns Columns to collapse into one condition column
#' @param condition_name Condition column name.
#' @returns A dataframe with additional columns condition_name, and condition_name + "_APPROX_EVENT_DAY"
multiple_columns_to_condition <- function(data,columns,condition_name="CONDITION"){
    condition_date = paste0(condition_name,"_APPROX_EVENT_DAY",sep="")
    columns_date = paste0(columns,"_APPROX_EVENT_DAY",sep="")
    dlist = list()
    for(i in 1:length(columns)){
        col = columns[[i]]
        d <- data %>% filter(get({{col}})==1)
        dlist[[i]] <- d
    }
    prelim_d <- do.call(rbind,dlist)
    prelim_d <- prelim_d %>% unique() %>% rowwise() %>% mutate(
        "{condition_name}" := max(c_across(all_of(columns)),na.rm=T),
        "{condition_date}" := min(c_across(all_of(columns_date)),na.rm=T)
    ) %>% ungroup()
    return(prelim_d)
}

#' Collapse multiple columns into a single event column, taking the first one as the event's date.
#'
#' @param data Dataframe
#' @param condition Condition column name. Will be used to filter the data, since the rowwise operations are slow.
#' @param columns Columns to collapse into one event column
#' @param event_name Event column name.
#' @returns A dataframe with additional columns event_name, and event_name + "_APPROX_EVENT_DAY"
multiple_columns_to_event <- function(data,condition,columns,event_name="EVENT"){
    event_date = paste0(event_name,"_APPROX_EVENT_DAY",sep="")
    columns_date = paste0(columns,"_APPROX_EVENT_DAY",sep="")
    prelim_d <- data %>% filter(get({{condition}})==1)
    prelim_d <- prelim_d %>% rowwise() %>% 
                            mutate("{event_name}" := max(c_across(all_of(columns)),na.rm=T),
                                   "{event_date}" := min(c_across(all_of(columns_date)),na.rm=T))%>%
                            ungroup()
    prelim_d[sapply(prelim_d,is.infinite)] <- NA
    return(prelim_d)
}