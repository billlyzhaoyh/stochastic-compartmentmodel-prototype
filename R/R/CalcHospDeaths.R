



##' 
##' Build a set of sampled hospitalizations, deaths, and recoveries 
##'  from the incident infection data from the simulation model.
##'  
##' @param data data.frame of t (time in days) and incidI (incident infections)
##' @param p_hosp probability of hospitalization, among infections
##' @param p_death probability of death, among infections (hospitalization is required for death)
##' @param time_hosp_pars parameters for time from onset to hospitalization distribution
##' @param time_death_pars parameters for time from hospitalization to death distribution
##' @param time_disch_pars parameters for time from hospitalization to discharge parameters
##' 
build_hospdeath_par_old <- function(data, p_hosp, p_death, p_hosp_type="gamma",
                            time_hosp_pars = c(1.23, 0.79), 
                            time_death_pars = c(log(11.25), log(1.15)), 
                            time_disch_pars = c(log(11.5), log(1.22)),
                            cores=4) {
    
    require(parallel)
    
    t_ <- 1:nrow(data)
    dates_ <- as.Date(data$date)
    uid <- data$sim_beta
    date_tmp <- seq(min(dates_), (max(dates_)+125), by="days")
    
    
    # Get counts of hospitalization, death, and recovery, by day
    I_ <- data$incid
    if (p_hosp_type=="single"){
        H_ <- rbinom(I_, I_, rep(p_hosp,  length(I_)))  # Add hosp
    } else if (p_hosp_type=="range"){
        H_ <- rbinom(I_, I_, runif(length(I_), p_hosp[1],  p_hosp[2]))  # Add hosp
    } else if (p_hosp_type=="gamma"){
        H_ <- rbinom(I_, I_, rgamma(length(I_), shape=p_hosp$shape, rate=p_hosp$rate))  # Add hosp
    }
    
    D_ <- rbinom(H_, H_, rep(p_death, length(H_))) # Add Death
    R_ <- H_ - D_                                  # Add hospitalized recoveries
    names(H_) <- uid
    names(D_) <- uid
    names(R_) <- uid
    
    
    # Time to hospitalization
    H_time_ <- floor(rlnorm(sum(H_), meanlog=time_hosp_pars[1], sdlog=time_hosp_pars[2]))
    H_date_ <- rep(dates_,H_) + H_time_
    names(H_time_) <- rep(names(H_), H_)
    names(H_date_) <- rep(names(H_), H_)
    
    # Time to death
    D_time_ <- floor(rlnorm(sum(D_), meanlog=time_death_pars[1], sdlog=time_death_pars[2]))
    D_date_ <- rep(dates_,D_) + D_time_
    D_start_ <- rep(dates_, D_)
    names(D_time_) <- rep(names(D_), D_)
    names(D_date_) <- rep(names(D_), D_)
    
    # Time to recovery
    R_time_ <- floor(rlnorm(sum(R_), meanlog=time_disch_pars[1], sdlog=time_disch_pars[2]))
    R_date_ <- rep(dates_,R_) + R_time_
    R_start_ <- rep(dates_,R_)
    names(R_time_) <- rep(names(R_), R_)
    names(R_date_) <- rep(names(R_), R_)
    
    
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    # seq_R_fn <- function(x=X, R_start_, R_date_){
    #     seq(as.Date(R_start_[x]), as.Date(R_date_[x]), "days")
    # }
    
    # Get current hospitalization days and accumulate them -- Recoveries
    clusterExport(cl=cl, varlist=c('R_start_', 'R_date_'), envir=environment())
    curr_hosp_date <- rev(as.Date(unlist(
        parSapply(cl, 1:sum(R_), function(x) seq(as.Date(R_start_[x]), as.Date(R_date_[x]), "days"))),
        origin = "1970-01-01"))
    names(curr_hosp_date) <- rep(names(R_time_), R_time_)
    
    # Get current hospitalization days and accumulate them -- Deaths
    clusterExport(cl=cl, varlist=c('D_start_', 'D_date_'), envir=environment())
    curr_hospD_date <- rev(as.Date(unlist(
        parSapply(cl, 1:sum(D_), function(x) seq(as.Date(D_start_[x]), as.Date(D_date_[x]), "days"))),
        origin = "1970-01-01"))
    names(curr_hospD_date) <- rep(names(D_time_), D_time_)
    
    stopCluster(cl)
    
    
    # combine them
    curr_hosp_date <- c(curr_hosp_date, curr_hospD_date)
    data_currhosp <- as.data.frame(table(curr_hosp_date, names(curr_hosp_date)))
    colnames(data_currhosp) <- c("t","sim_beta","hosp_curr")
    
    
    data_H <- as.data.frame(table(H_date_, names(H_date_)))
    data_D <- as.data.frame(table(D_date_, names(D_date_)))
    colnames(data_H) <- c("t", "sim_beta", "incidH")
    colnames(data_D) <- c("t", "sim_beta", "incidD")
    
    data_D$t <- as.Date(data_D$t)
    data_H$t <- as.Date(data_H$t)
    data_currhosp$t <- as.Date(data_currhosp$t)
    
    
    res <- full_join(data.frame(t=as.Date(date_tmp)), data_H, by=c("t"="t"))
    res <- full_join(res, data_D, by=c("t"="t", "sim_beta"="sim_beta"))
    res <- full_join(res, data_currhosp, by=c("t"="t", "sim_beta"="sim_beta"))
    
    res <- res %>% mutate(incidH = ifelse(is.na(incidH), 0, incidH),
                          incidD = ifelse(is.na(incidD), 0, incidD),
                          hosp_curr = ifelse(is.na(hosp_curr), 0, hosp_curr))
    
    res <- res %>% filter(!is.na(sim_beta)) %>%
        arrange(sim_beta, t)
    
    return(res)
}










# RUN THE CALCULATION -----------------------------------------------------

# *** will have to adjust to model output and time ****
# *** add detection and confirmation later
# 
# test_data <- data.frame(t = 0:8, S = c(1000, 999, 999, 990, 980, 960, 800, 700, 500), incidI = c(1, 1, 0, 10, 10, 20, 160, 100, 200))
# 
# 
# data <- test_data %>% select(t, incidI)
# res_data <- build_hospdeath(data, p_hosp=p_hosp[3], p_death=p_death[3], time_hosp_pars, time_death_pars, time_disch_pars)


create_delay_frame <- function(X, p_X, data_, X_pars, varname, p_X_type="single") {
    
    if (p_X_type=="single"){
        X_ <- rbinom(length(data_[[X]]),data_[[X]],p_X)
    } else if (p_X_type=="range"){
        X_ <- rbinom(data_[[X]], data_[[X]], runif(length(data_[[X]]), p_X[1],  p_X[2]))
    } else if (p_X_type=="gamma"){
        X_ <- rbinom(data_[[X]], data_[[X]], rgamma(length(data_[[X]]), shape=p_X$shape, rate=p_X$rate))
    }
    
    #X_ <- rbinom(length(data_[[X]]),data_[[X]],p_X)
    data_X <- data.frame(time=data_$time,  uid=data_$uid, count=X_)
    X_delay_ <- round(exp(X_pars[1] + X_pars[2]^2 / 2))
    
    X_time_ <- rep(as.Date(data_X$time), data_X$count) + X_delay_
    names(X_time_) <- rep(data_$uid, data_X$count)
    
    data_X <- data.frame(time=X_time_, uid=names(X_time_))
    data_X <- data.frame(setDT(data_X)[, .N, by = .(time, uid)])
    colnames(data_X) <- c("time","uid",paste0("incid",varname))
    return(data_X)
}




##' 
##' Build a set of sampled hospitalizations, deaths, and recoveries 
##'  from the incident infection data from the simulation model.
##'  
##' @param data data.frame of t (time in days) and incidI (incident infections)
##' @param p_hosp probability of hospitalization, among infections
##' @param p_death probability of death, among infections (hospitalization is required for death)
##' @param time_hosp_pars parameters for time from onset to hospitalization distribution
##' @param time_death_pars parameters for time from hospitalization to death distribution
##' @param time_disch_pars parameters for time from hospitalization to discharge parameters
##' 
build_hospdeath_par <- function(data, p_hosp, p_death, p_vent, p_ICU, p_hosp_type="gamma",
                                time_hosp_pars = c(1.23, 0.79), 
                                time_ICU_pars = c(log(10.5), log((10.5-7)/1.35)),
                                time_vent_pars = c(log(10.5), log((10.5-8)/1.35)),
                                time_death_pars = c(log(11.25), log(1.15)), 
                                time_disch_pars = c(log(11.5), log(1.22)),
                                time_ICUdur_pars = c(log(16), log(2.96)),
                                end_date = "2022-01-01",
                                length_geoid = 5,
                                cores=4, 
                                run_parallel=FALSE){
    
    library(data.table)
    library(doParallel)
    
    dat_final <- list()
    
    
    n_sim <- length(unique(data$sim_num))
    print(paste("Creating cluster with",cores,"cores"))
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    
    print(paste("Running over",n_sim,"simulations"))
    dat_final <- foreach(s=seq_len(n_sim), .packages=c("dplyr","readr","data.table","tidyr"),
                         .export = c("create_delay_frame")) %dopar% {
                             print(s)
                             
                             dat_ <- data %>% filter(sim_num==s)
                             # county_dat <- read.csv("data/west-coast-AZ-NV/geodata.csv")
                             # county_dat$geoid <- as.character(county_dat$geoid)
                             # county_dat$new_pop <- county_dat$pop2010
                             # county_dat <- make_metrop_labels(county_dat)
                             # dat_ <- load_scenario_sim(data_filename,s,keep_compartments = c("diffI","cumI")) %>%
                             #     filter(geoid %in% county_dat$geoid[county_dat$stateUSPS=="CA"], time<=end_date, comp == "diffI", N > 0) %>%
                             #     mutate(hosp_curr = 0, icu_curr = 0, vent_curr = 0, uid = paste0(geoid, "-",sim_num)) %>%
                             #     rename(incidI = N)
                             dates_ <- as.Date(dat_$time)
                             
                             # Add time things
                             dat_H <- create_delay_frame('incidI',p_hosp,dat_,time_hosp_pars,"H", p_X_type = "gamma")
                             data_ICU <- create_delay_frame('incidH',p_ICU,dat_H,time_ICU_pars,"ICU")
                             data_Vent <- create_delay_frame('incidICU',p_vent,data_ICU,time_vent_pars,"Vent")
                             data_D <- create_delay_frame('incidH',p_death,dat_H,time_death_pars,"D")
                             
                             
                             R_delay_ <- round(exp(time_disch_pars[1]))
                             ICU_dur_ <- round(exp(time_ICUdur_pars[1]))
                             
                             
                             
                             # Using `merge` instead     
                             res <- merge(dat_H %>% mutate(uid = as.character(uid)), 
                                          data_ICU %>% mutate(uid = as.character(uid)), all=TRUE)
                             res <- merge(res, data_Vent %>% mutate(uid = as.character(uid)), all=TRUE)
                             res <- merge(res, data_D %>% mutate(uid = as.character(uid)), all=TRUE)
                             res <- merge(dat_ %>% mutate(uid = as.character(uid)), 
                                          res %>% mutate(uid = as.character(uid)), all=TRUE)
                             
                             res <- res %>% 
                                 replace_na(
                                     list(incidI = 0,
                                          incidH = 0,
                                          incidICU = 0,
                                          incidVent = 0,
                                          incidD = 0,
                                          vent_curr = 0,
                                          hosp_curr = 0))
                             
                             # get sim nums
                             res <- res %>% select(-geoid, -sim_num) %>%
                                 separate(uid, c("geoid", "sim_num"), sep="-", remove=FALSE)
                             
                             res <- res %>% mutate(date_inds = as.integer(time - min(time) + 1))
                             n_sim <- length(unique(res$sim_num))
                             
                             
                             
                             res$sim_num_good <- as.numeric(res$sim_num) 
                             res$sim_num_good <- res$sim_num_good - min(res$sim_num_good) +1
                             
                             res$geo_ind <- as.numeric(as.factor(res$geoid))
                             inhosp <- matrix(0, nrow=max(res$date_inds),ncol=max(res$geo_ind))
                             inicu <- inhosp
                             len<-max(res$date_inds)
                             for (i in 1:nrow(res)) {
                                 inhosp[res$date_inds[i]:min((res$date_inds[i]+R_delay_-1),len), res$geo_ind[i]] <- 
                                     inhosp[res$date_inds[i]:min((res$date_inds[i]+R_delay_-1),len), res$geo_ind[i]] + res$incidH[i]
                                 
                                 inicu[res$date_inds[i]:min((res$date_inds[i]+ICU_dur_-1),len), res$sim_num_good[i]] <- 
                                     inicu[res$date_inds[i]:min((res$date_inds[i]+ICU_dur_-1),len), res$sim_num_good[i]] + res$incidICU[i]
                                 
                             }
                             
                             
                             
                             for (x in 1:nrow(res)){
                                 
                                 res$hosp_curr[x] <- inhosp[res$date_inds[x], res$geo_ind[x]]
                                 res$icu_curr[x] <- inicu[res$date_inds[x], res$geo_ind[x]]
                                 #res$hosp_curr <- res$hosp_curr + res$date_inds %in% (res$date_inds[x] + 0:R_delay_)*res$incidH[x]
                                 #res$icu_curr <- res$icu_curr + res$date_inds %in% (res$date_inds[x] + 0:ICU_dur_)*res$incidICU[x]
                                 #res$vent_curr <- res$vent_curr + res$date_inds %in% (res$date_inds[x] + 0:Vent_dur_)*res$incidVent[x]
                             }
                             
                             dat_final[[s]] <- res
                             
                         }
    
    dat_final <- data.table::rbindlist(dat_final)
    
    print(paste("Parallel portion finished"))
    stopCluster(cl)
    
    return(dat_final)
}







