library(tidyverse)

##' estimates probability of being a severe case by age from Shenzen data
##' results from this STAN model saved for easy use later on
get_severe_age_shenzhen <- function( ){

    # Raw data not able to be shared
    return("raw data not available to the public.")
    

    dat <- full_join(sym_dat, fev_dat)
    dat <- dat %>% as.data.frame() %>% mutate(tot = yes+no)
    dat <- dat %>% mutate(not_severe = tot-severe)
    dat <- dat %>% mutate(p_severe = severe/tot)

    # Use stan

    library(rstanarm)

    t_prior <- student_t(df = 7, location = 0, scale = 2.5, autoscale = FALSE)
    fit1 <- stan_glm(cbind(severe, tot-severe) ~ age_cat, data = dat,
                     family = binomial(link = "logit"),
                     prior = t_prior, prior_intercept = t_prior,
                     cores = 4, seed = 12345)

    PPD <- posterior_predict(fit1)
    prob <- PPD
    for(i in 1:nrow(PPD)){
        prob[i,] <- PPD[i,] / dat$tot
    }

    write_csv(prob %>% as.data.frame(), "data/severe_age_prob.csv")
    return(prob)

}



##'
##' Get population distribution and aggregate it to 10 year age groups
##' - this is set up to use population estimates from the World Populaiton Prospects data
##'
##' @param country country of interest
##'
get_age_pop <- function(country){

    require(stringi)
    #require(globaltoolbox)

    pop_data <- read_csv("data/WPP2019_POP.csv")
    pop_data <- pop_data %>%
        mutate(country_clean = stringi::stri_trans_general(location, "Latin-ASCII")) %>%
        filter(tolower(country_clean) == tolower(country)) %>% filter(year==max(year))

    # print for a double check
    print(pop_data$location)
    pop_data <- pop_data[,-(1:4)] %>% select(-country_clean)
    dat <- as.numeric(pop_data)
    names(dat) <- colnames(pop_data)
    return(dat)
}





##'
##' Get population distribution and aggregate it to 10 year age groups
##'  - this is set up to use population estimates from the World Populaiton Prospects data
##'
##' @param country country of interest
##'
get_p_severe <- function(country="China"){

    # Load prob(severe | age) from shenzhen
    prob <- read_csv("data/severe_age_prob.csv")


    #  population by age
    nage_ <- get_age_pop(country) * 1000
    nage_[8] <- sum(nage_[8:11])
    nage_ <- nage_[1:8]
    pr_age10_ <- nage_ / sum(nage_)

    p_severe_tmp <- prob
    for(i in 1:nrow(prob)){
        p_severe_tmp[i,] <- prob[i,] * pr_age10_
    }
    p_severe_ <- rowSums(p_severe_tmp)

    fit_ <- fitdistrplus::fitdist(p_severe_, "gamma", "mle")


    p_severe_ <- list(ests = p_severe_,
                      mean=mean(p_severe_),
                      ll=quantile(p_severe_, .025),
                      ul=quantile(p_severe_, .975),
                      q25=quantile(p_severe_, .25),
                      q75=quantile(p_severe_, .75),
                      shape = coef(fit_)["shape"],
                      rate = coef(fit_)["rate"])

    return(p_severe_)
}



##' estimates proportion of severe cases adjusted for population structure
##'
##' @param pr_age10 proportion of population in 10 year age bins
##'
get_p_severe_pop <- function(pr_age10){

    # Load prob(severe | age) from shenzhen
    prob <- read_csv("data/severe_age_prob.csv")

    #  sum all proportion of age old than 70
    pr_age10[8] <- sum(pr_age10[8:length(pr_age10)])

    p_severe_tmp <- prob
    for(i in 1:nrow(prob)){
        p_severe_tmp[i,] <- prob[i,] * pr_age10
    }
    p_severe_ <- rowSums(p_severe_tmp)

    fit_ <- fitdistrplus::fitdist(p_severe_, "gamma", "mle")


    p_severe_ <- list(ests = p_severe_,
                      mean=mean(p_severe_),
                      ll=quantile(p_severe_, .025),
                      ul=quantile(p_severe_, .975),
                      q25=quantile(p_severe_, .25),
                      q75=quantile(p_severe_, .75),
                      shape = coef(fit_)["shape"],
                      rate = coef(fit_)["rate"])

    return(p_severe_)
}
