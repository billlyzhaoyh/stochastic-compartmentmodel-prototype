
library(doParallel)

##' Model Functions
##' 
##' 
run_seir_model <- function(n.sims=20, N=600000, betas=betas_low, sigma=sigma, gamma=gamma, 
                           seeds=1, start_date_="2020-03-15", step.size=.25, run_time=360){
    
    sim_res <- list()
    
    for (b in 1:length(betas)){
        
        sim_res_sim <- list()
        
        for (i in 1:n.sims){
            tmp <- runStochasticSEIR(beta = betas[b],
                                    sigma = sigma,
                                    gamma = gamma,
                                    initial.state = c(S=N-seeds, E=seeds,I1=0,I2=0,I3=0,R=0,I_trk=0,E_trk=0),
                                    step.size = step.size,
                                    freq.dependent = TRUE,
                                    final.only = FALSE,
                                    max.time = run_time)
            
            # get epi curve by aggregtating 
            new.cases <- diff(tmp[,"I_trk"])
            new.exposed <- diff(tmp[,"E_trk"])
            
            days <- floor(tmp[,"time"])
            days_m1 <- days[-1]
            epi.curve <- tapply(new.cases,days_m1,sum) # sum up time steps by day
            inf.curve <- tapply(new.exposed,days_m1,sum) # sum up time steps by day
            
            sim_res_sim[[i]] <- data.frame(sim=i, beta_ind=b, beta=betas[b], 
                                           day=unique(days), seeds=seeds,
                                           incid=epi.curve, 
                                           exposed=inf.curve,
                                           cum_exposed=cumsum(inf.curve),
                                           cum_cases=cumsum(epi.curve))
        }
        
        sim_res[[b]] <- data.table::rbindlist(sim_res_sim)
        
    }
    
    sim_res <- data.table::rbindlist(sim_res)
    
    dates <- seq.Date(as.Date(start_date_), as.Date(start_date_)+run_time-1, by=1)
    sim_res$date <- start_date_ + sim_res$day
    sim_res <- sim_res %>% arrange(sim, beta) %>% 
        mutate(sim_beta = paste0(sim, "-", beta_ind)) %>% 
        mutate(sim_id = as.integer(as.factor(sim_beta)))
    
    return(sim_res)
}








run_seir_model_parallel <- function(project, n.sims=20, N=600000, betas=betas_low, R0_val="low", sigma=sigma, gamma=gamma, 
                           seeds=1, start_date_="2020-03-15", step.size=.25, run_time=360, cores=4){
    
    sim_res <- list()
    
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    for (b in 1:length(betas)){
        
        sim_res_sim <- foreach(i=1:n.sims, .export = c("runStochasticSEIR", "stochastic.dx.dt")) %dopar% {
            
            tmp <- runStochasticSEIR(beta = betas[b],
                                     sigma = sigma,
                                     gamma = gamma,
                                     initial.state = c(S=N-seeds, E=seeds,I1=0,I2=0,I3=0,R=0,I_trk=0,E_trk=0),
                                     step.size = step.size,
                                     freq.dependent = TRUE,
                                     final.only = FALSE,
                                     max.time = run_time)
            
            # get epi curve by aggregtating 
            new.cases <- diff(tmp[,"I_trk"])
            new.exposed <- diff(tmp[,"E_trk"])
            
            days <- floor(tmp[,"time"])
            days_m1 <- days[-1]
            epi.curve <- tapply(new.cases,days_m1,sum) # sum up time steps by day
            inf.curve <- tapply(new.exposed,days_m1,sum) # sum up time steps by day
            
            data.frame(sim=i, beta_ind=b, beta=betas[b], 
                        day=unique(days), seeds=seeds,
                        incid=epi.curve, 
                        exposed=inf.curve,
                        cum_exposed=cumsum(inf.curve),
                        cum_cases=cumsum(epi.curve))
        }
        
        sim_res[[b]] <- data.table::rbindlist(sim_res_sim)
        
    }
    
    stopCluster(cl)
    
    sim_res <- data.table::rbindlist(sim_res)
    
    dates <- seq.Date(as.Date(start_date_), as.Date(start_date_)+run_time-1, by=1)
    sim_res$date <- as.Date(start_date_) + sim_res$day
    sim_res <- sim_res %>% arrange(sim, beta) %>% 
        mutate(sim_beta = paste0(sim, "-", beta_ind)) %>% 
        mutate(sim_id = as.integer(as.factor(sim_beta)))
    
    
    dir.create(file.path("results",project), recursive = TRUE, showWarnings = FALSE)
    write_csv(sim_res, file.path("results", project, paste0("sim_res_nohosp_", R0_val,"-", 
                                         month.abb[lubridate::month(Sys.Date())],lubridate::day(Sys.Date()), ".csv")))
    return(sim_res)

}




