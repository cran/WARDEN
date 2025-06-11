## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE
)

## ----setup--------------------------------------------------------------------
library(WARDEN)

library(dplyr)
library(ggplot2)
library(kableExtra)
library(purrr)


## ----main_opt, results='hide', message=FALSE----------------------------------
options(scipen = 999)
options(digits=3)
options(tibble.print_max = 50)

## ----input_delayed------------------------------------------------------------
#We don't need to use sensitivity_inputs here, so we don't add that object

#Put objects here that do not change on any patient or intervention loop
#We use add_item2 and add_item to showcase how the user can implement the inputs (either works, add_item2 is just faster)
common_all_inputs <-add_item2(input = {
                      util.sick   <- 0.8
                      util.sicker <- 0.5
                      cost.sick   <- 3000
                      cost.sicker <- 7000
                      cost.int    <- 1000
                      coef_noint  <- log(0.2)
                      HR_int      <- 0.8
                      drc         <- 0.035 #different values than what's assumed by default
                      drq         <- 0.035
                      random_seed_sicker_i <- sample.int(100000,npats,replace = FALSE)
})  #to be used as seeds to draw the time to event for sicker, to ensure same luck for the same patient independently of the arm


#Put objects here that do not change as we loop through treatments for a patient
common_pt_inputs <- add_item(death= max(0.0000001,rnorm(n=1, mean=12, sd=3))) 

#Put objects here that change as we loop through treatments for each patient (e.g. events can affect fl.tx, but events do not affect nat.os.s)
unique_pt_inputs <- add_item(fl.sick = 1,
                             q_default = util.sick,
                             c_default = cost.sick + if(arm=="int"){cost.int}else{0}) 


## ----model_evts---------------------------------------------------------------
init_event_list <- 
  add_tte(arm=c("noint","int"), evts = c("sick","sicker","death") ,input={
    sick <- 0
    sicker <- draw_tte(1,dist="exp", coef1=coef_noint, beta_tx = ifelse(arm=="int",HR_int,1), seed = random_seed_sicker_i[i]) #this way the value would be the same if it wasn't for the HR, effectively "cloning" patients luck
    
  })


## ----model_reaction-----------------------------------------------------------
evt_react_list <-
  add_reactevt(name_evt = "sick",
               input = {}) %>%
  add_reactevt(name_evt = "sicker",
               input = {
                 q_default <- util.sicker
                 c_default <- cost.sicker + if(arm=="int"){cost.int}else{0}
                 fl.sick   <- 0 
               }) %>%
  add_reactevt(name_evt = "death",
               input = {
                 q_default <- 0
                 c_default <- 0
                 curtime   <- Inf
               }) 


# Below how it would be set up if using `accum_backwards = TRUE` in `run_sim()` (and will give equal final results)
# Note that we set the value applied in the reaction right up to the event, changing the interpretation of the reaction
# It is also a slower method than the standard approach
# In this case we use `modify_item_seq` in the reaction for death, but given that the items are not interacting with each other,
# we could use `modify_item` instead and make the computation faster
#
# evt_react_list <-
#     add_reactevt(name_evt = "sick",
#                  input = {}) %>%
#     add_reactevt(name_evt = "sicker",
#                  input = {
#                      modify_item(list(
#                      q_default = util.sick,
#                      c_default = cost.sick + if(arm=="int"){cost.int}else{0},
#                      fl.sick = 0
#                      ))
#                  }) %>%
#     add_reactevt(name_evt = "death",
#                  input = {
#                      c_dis <- if(fl.sick==1){cost.sick}else{cost.sicker} 
#                      modify_item_seq(list(
#                      q_default = if(fl.sick==1){util.sick}else{util.sicker},
#                      c_default = c_dis + if(arm=="int"){cost.int}else{0},
#                      curtime = Inf
#                      ))
#                  }) 


## ----model_relatioships-------------------------------------------------------

df_interactions <- extract_from_reactions(evt_react_list)

kable(df_interactions)


## ----utilities----------------------------------------------------------------

util_ongoing <- "q_default"


## ----costs--------------------------------------------------------------------

cost_ongoing <- "c_default"


## ----model_run----------------------------------------------------------------
#Logic is: per patient, per intervention, per event, react to that event.
results <- run_sim(  
  npats=1000,                               # number of patients to be simulated
  n_sim=1,                                  # number of simulations to run
  psa_bool = FALSE,                         # use PSA or not. If n_sim > 1 and psa_bool = FALSE, then difference in outcomes is due to sampling (number of pats simulated)  
  arm_list = c("int", "noint"),             # intervention list
  common_all_inputs = common_all_inputs,    # inputs common that do not change within a simulation
  common_pt_inputs = common_pt_inputs,      # inputs that change within a simulation but are not affected by the intervention
  unique_pt_inputs = unique_pt_inputs,      # inputs that change within a simulation between interventions
  init_event_list = init_event_list,        # initial event list
  evt_react_list = evt_react_list,          # reaction of events
  util_ongoing_list = util_ongoing,
  cost_ongoing_list = cost_ongoing,
  ipd = 1
)

## ----post-processing_summary--------------------------------------------------


summary_results_det(results[[1]][[1]]) #print first simulation

summary_results_sim(results[[1]])

summary_results_sens(results)

psa_ipd <- bind_rows(map(results[[1]], "merged_df")) 

psa_ipd[1:10,] %>%
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))



## ----post-processing_analysis,echo=FALSE, message=FALSE-----------------------

psa_ipd %>% group_by(arm,evtname) %>% summarise(n=n()) %>% arrange(arm,-n)%>%
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))


## ----post-processing_plots1, fig.width=10, fig.height=8-----------------------

data_plot <- results[[1]][[1]]$merged_df %>%
  filter(evtname != "sick") %>%
  group_by(arm,evtname,simulation) %>%
  mutate(median = median(evttime)) %>%
  ungroup()

ggplot(data_plot) +
  geom_density(aes(fill = arm, x = evttime),
               alpha = 0.7) +
  geom_vline(aes(xintercept=median,col=arm)) +
  facet_wrap( ~ evtname, scales = "free") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_bw()


## ----post-processing_plots, fig.width=10, fig.height=8, message=FALSE---------

data_qaly_cost<- psa_ipd[,.SD[1],by=.(pat_id,arm,simulation)][,.(arm,qaly=total_qalys,cost=total_costs,pat_id,simulation)]
data_qaly_cost[,ps_id:=paste(pat_id,simulation,sep="_")]


mean_data_qaly_cost <- data_qaly_cost %>% group_by(arm) %>% summarise(across(where(is.numeric),mean))

ggplot(data_qaly_cost,aes(x=qaly, y = cost, col = arm)) + 
  geom_point(alpha=0.15,shape = 21) +
  geom_point(data=mean_data_qaly_cost, aes(x=qaly, y = cost, fill = arm), shape = 21,col="black",size=3) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = .5))


## ----dsa_inputs---------------------------------------------------------------
#Load some data
list_par <- list(parameter_name = list("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","HR_int"),
                              base_value = list(0.8,0.5,3000,7000,1000,log(0.2),0.8),
                              DSA_min = list(0.6,0.3,1000,5000,800,log(0.1),0.5),
                              DSA_max = list(0.9,0.7,5000,9000,2000,log(0.4),0.9),
                              PSA_dist = list("rnorm","rbeta_mse","rgamma_mse","rgamma_mse","rgamma_mse","rnorm","rlnorm"),
                              a=list(0.8,0.5,3000,7000,1000,log(0.2),log(0.8)),
                              b=lapply(list(0.8,0.5,3000,7000,1000,log(0.2),log(0.8)), function(x) abs(x/5)),
                              scenario_1=list(0.6,0.3,1000,5000,800,log(0.1),0.5),
                              scenario_2=list(0.9,0.7,5000,9000,2000,log(0.4),0.9)
                              )

sensitivity_inputs <-add_item(
            indicators = if(sensitivity_bool){ create_indicators(sens,n_sensitivity*length(sensitivity_names),rep(1,length(list_par[[1]])))}else{
                                rep(1,length(list_par[[1]]))} #vector of indicators, value 0 everywhere except at sens, where it takes value 1 (for dsa_min and dsa_max, if not sensitivity analysis, then we activate all of them, i.e., in a PSA)
                              )

common_all_inputs <-add_item() %>% 
  add_item(
            pick_val_v(base        = list_par[["base_value"]],
                       psa         = pick_psa(list_par[["PSA_dist"]],rep(1,length(list_par[["PSA_dist"]])),list_par[["a"]],list_par[["b"]]),
                       sens        = list_par[[sens_name_used]],
                       psa_ind     = psa_bool,
                       sens_ind    = sensitivity_bool,
                       indicator   = indicators,
                       names_out   = list_par[["parameter_name"]]
                       )
            ) %>%
  add_item(random_seed_sicker_i = sample(1:1000,1000,replace = FALSE)) #we don't add this variable ot the sensitivity analysis


## ----run_dsa------------------------------------------------------------------
results <- run_sim(  
  npats=100,                               # number of patients to be simulated
  n_sim=1,                                  # number of simulations to run
  psa_bool = FALSE,                         # use PSA or not. If n_sim > 1 and psa_bool = FALSE, then difference in outcomes is due to sampling (number of pats simulated)  
  arm_list = c("int", "noint"),             # intervention list
  common_all_inputs = common_all_inputs,    # inputs common that do not change within a simulation
  common_pt_inputs = common_pt_inputs,      # inputs that change within a simulation but are not affected by the intervention
  unique_pt_inputs = unique_pt_inputs,      # inputs that change within a simulation between interventions
  init_event_list = init_event_list,        # initial event list
  evt_react_list = evt_react_list,          # reaction of events
  util_ongoing_list = util_ongoing,
  cost_ongoing_list = cost_ongoing,
  sensitivity_inputs = sensitivity_inputs,
  sensitivity_names = c("DSA_min","DSA_max"),
  sensitivity_bool = TRUE,
  n_sensitivity = length(list_par[[1]]),
  input_out = unlist(list_par[["parameter_name"]])
)


## ----dsa_check----------------------------------------------------------------

data_sensitivity <- bind_rows(map_depth(results,2, "merged_df"))

#Check mean value across iterations as PSA is off
data_sensitivity %>% group_by(sensitivity) %>% summarise_at(c("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","HR_int"),mean)


## ----run_dsa_psa--------------------------------------------------------------
results <- run_sim(  
  npats=100,                               
  n_sim=6,                                  
  psa_bool = TRUE,                         
  arm_list = c("int", "noint"),             
  common_all_inputs = common_all_inputs,    
  common_pt_inputs = common_pt_inputs,      
  unique_pt_inputs = unique_pt_inputs,      
  init_event_list = init_event_list,        
  evt_react_list = evt_react_list,          
  util_ongoing_list = util_ongoing,
  cost_ongoing_list = cost_ongoing,
  sensitivity_inputs = sensitivity_inputs,
  sensitivity_names = c("DSA_min","DSA_max"),
  sensitivity_bool = TRUE,
  n_sensitivity = length(list_par[[1]]),
  input_out = unlist(list_par[["parameter_name"]])
)


## ----dsa_check_psa------------------------------------------------------------

data_sensitivity <- bind_rows(map_depth(results,2, "merged_df"))

#Check mean value across iterations as PSA is off
data_sensitivity %>% group_by(sensitivity) %>% summarise_at(c("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","HR_int"),mean)


## ----run_psa------------------------------------------------------------------
results <- run_sim(  
  npats=100,                               
  n_sim=10,                                  
  psa_bool = TRUE,                         
  arm_list = c("int", "noint"),             
  common_all_inputs = common_all_inputs,    
  common_pt_inputs = common_pt_inputs,      
  unique_pt_inputs = unique_pt_inputs,      
  init_event_list = init_event_list,        
  evt_react_list = evt_react_list,          
  util_ongoing_list = util_ongoing,
  cost_ongoing_list = cost_ongoing,
  sensitivity_inputs = sensitivity_inputs,
  sensitivity_bool = FALSE,
  n_sensitivity = 1,
  input_out = unlist(list_par[["parameter_name"]])
)


## ----check_psa----------------------------------------------------------------

data_simulation <- bind_rows(map_depth(results,2, "merged_df"))

#Check mean value across iterations as PSA is off
data_simulation %>% group_by(simulation) %>% summarise_at(c("util.sick","util.sicker","cost.sick","cost.sicker","cost.int","coef_noint","HR_int"),mean)


