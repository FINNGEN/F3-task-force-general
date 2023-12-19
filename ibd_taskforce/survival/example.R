source("./functions.R")

#input your favourite datapath
datapath = ""
data <- read_tsv(datapath,guess_max=1000000)

condition_col  = ""
event_col = ""
censor_col = ""
start_date = "1995-01-01"
end_date = "2015-01-01"
followup_period = 5
additional_columns=c()


survival_data <- survival_analysis_helper(data,
                                        condition_col,
                                        event_col,
                                        start_date,
                                        end_date,
                                        followup_period,
                                        censor_col,
                                        additional_columns,
                                        )
survival_fit <- coxph(Surv(time,status)~date_modifier,data=survival_data)
summary(survival_fit)
# for plotting purposes, the stratified factor can be useful
survival_fit2 <- survfit(Surv(time,status)~date_interval,data=survival_data)
ggsurvplot(survival_fit2,conf.int=T,censor=F)

# If you want to combine multiple columns to be either the condition or event column, that can be done
# These helpers take the first event (as in first occurred event, in dates) as the event day
multiple_condition_columns <- c("col1","col2")
condition_name <- "CONDITION"
new_data <- multiple_columns_to_condition(data,multiple_condition_columns,condition_name)
# the columns are then available as CONDITION and CONDITION_APPROX_EVENT_DAY
survival_data <- survival_data <- survival_analysis_helper(new_data,
                                        condition_name,
                                        event_col,
                                        start_date,
                                        end_date,
                                        followup_period,
                                        censor_col,
                                        additional_columns,
                                        )

# Or to combine multiple columns to one event column, that can be done
# Again, the event date will be the first of the events that happened.
multiple_event_columns = c("col1","col2")
new_event_name <- "EVENT"
new_data <- multiple_columns_to_event(data,condition,multiple_event_columns,new_event_name)
survival_data <- survival_data <- survival_analysis_helper(new_data,
                                        condition_name,
                                        new_event_name,
                                        start_date,
                                        end_date,
                                        followup_period,
                                        censor_col,
                                        additional_columns,
                                        )