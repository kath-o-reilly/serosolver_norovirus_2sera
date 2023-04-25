######################################################
## SIMULATION RECOVERY TEST -- NOROVIRUS DATA WITH BINDING AVIDITY
## Author: James Hay
## Date: 13 April 2023
## Summary: simulates some serosurvey data with two "observation types" and fits serosolver

library(ggplot2)
library(coda)
library(plyr)
library(reshape2)
library(data.table)
library(foreach)
library(doParallel)
library(parallel)
library(dplyr)
library(readr)
library(viridis)
library(ggpubr)
library(tidyverse)

serosolver_wd <- "~/GitHub/serosolver/"
devtools::load_all(serosolver_wd)
##library(serosolver) 

packageDescription("serosolver") # check it is the right branch

run_name <- "sim_noro"
main_wd <- "~/GitHub/serosolver_norovirus_2sera/"
chain_wd <- paste0(main_wd,"/chains_avidity/",run_name)
save_wd <- paste0(main_wd,"/figures/chain_plots_avidity/")

if(!dir.exists(save_wd)) dir.create(save_wd,recursive = TRUE)
if(!dir.exists(chain_wd)) dir.create(chain_wd,recursive = TRUE)

buckets <- 1 ## Ignore
prior_version <- 2 ## Which prior on the infection histories to use? Prior version 2 is generally preferred
n_chains <- 5 ## Number of MCMC chains to run

rerun <- TRUE ## Set to FALSE if you just want to load in previously run chains

setwd(main_wd)
print(paste0("In directory: ", main_wd))
print(paste0("Saving to: ", save_wd))

set.seed(1)

## Run multiple chains in parallel
cl <- makeCluster(n_chains)
registerDoParallel(cl)

## MCMC settings, not super important but can be tweaked
mcmc_pars <- c("save_block"=100,"thin"=10,"thin_hist"=50,
               "iterations"=200000,
               "adaptive_period"=100000,
               "burnin"=0,"switch_sample"=2,"hist_switch_prob"=0.05,
               "year_swap_propn"=0.8,"swap_propn"=0.5,
               "inf_propn"=0.5,"hist_sample_prob"=1,"move_size"=3, "hist_opt"=0,
               "popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,propose_from_prior=TRUE)

## Simulation parameters
n_obs_types <- 2
obs_type_dists <- c(1,2) ## Vector of observation models for each observation type -- use 1 for discretized normal (as in previous version) and 2 for continuous, truncated normal
n_indivs <- 250
n_groups <- 1 ## Leave as 1, can in theory make multiple groups with distinct FOI parameters, but needs extra setup
n_samps <- 2 ## Number of samples per person
repeats <- 1 ## Number of repeat measurements per variant/sample combination
samp_min <- 2009 ## First sample year
samp_max <- 2012 ## Final sample year
year_min <- 2000 ## First year of possible circulation (ie. time 0 of the simulation)
year_max <- 2012 ## Final year of possible circulation
age_min <- 2 ## Age minimum and maximum in years, simulated from a uniform distribution
age_max <- 10

## Viruses and times for samples
sampled_viruses <- c(2000,2002,2006,2009,2012)
sampling_times <- seq(samp_min, samp_max, by=1)

## Create a fake antigenic map -- can put in what you like here
antigenic_coords <- data.frame(Strain=c(2000,2002,2006,2009,2012),X=c(0,0.5,3,3.5,4),Y=c(0,2,1,3,4))
antigenic_map <- generate_antigenic_map_flexible(antigenic_coords,
                                                 year_max=2013,year_min=2000,spar = 0.001)
ggplot(antigenic_map) + geom_line(aes(x=x_coord,y=y_coord)) +
    geom_point(data=antigenic_coords,aes(x=X,y=Y))

strain_isolation_times <- antigenic_map$inf_times
n_times <- length(strain_isolation_times)

## Set up parameter table
par_tab <- read.csv("par_tab_base_new.csv",stringsAsFactors=FALSE)
par_tab <- par_tab[par_tab$names != "phi",]
par_tab[par_tab$names %in% c("alpha","beta"),c("values")] <- c(1/3,1/3) ## Can also try c(1,1), or something informative. Just the parameters of a beta distribution which acts as the prior on the per-time attack rate.

## Just some setup of the parameter table to get parameters vaguely like you showed me
par_tab$fixed <- 1
par_tab[par_tab$names %in% c("mu","tau","sigma1","obs_sd","mu_short","wane","sigma2"),"fixed"] <- 0

## Create some random measurement offsets -- these add a systematic shift to each observed variant. Only those sampled variant years have offsets, the rest are 0
## These are unknown parameters to be estimated, assuming the offsets are drawn from ~ norm(0, 1)
par_tab_rhos <- data.frame(names="rho",values=rep(0,length(sampled_viruses)),fixed=0,steps=0.1,
                           lower_bound=-3,upper_bound=3,lower_start=-1,
                           upper_start=1,type=3)

par_tab_rhos$values <- rnorm(length(sampled_viruses), 1)
measurement_indices <- seq_along(sampled_viruses)

measurement_indices <- data.frame(virus = sampled_viruses, 
                                  obs_type = rep(1:n_obs_types, each=length(sampled_viruses)),
                                  rho_index=1:(length(sampled_viruses)*n_obs_types))

par_tab <- bind_rows(par_tab, par_tab_rhos)

## Extend parameter table for each aditional observation type
par_tab$obs_type <- 1
antigenic_map_tmp <- antigenic_map
antigenic_map$obs_type <- 1
par_tab_tmp <- par_tab
if(n_obs_types > 1){
    for(i in 2:n_obs_types){
        par_tab_tmp2 <- par_tab_tmp
        antigenic_map_tmp2 <- antigenic_map_tmp
        antigenic_map_tmp2$obs_type <- i 
        par_tab_tmp2$obs_type <- i
        par_tab <- bind_rows(par_tab, par_tab_tmp2 %>% filter(!(names %in% c("alpha","beta"))))
        antigenic_map <- bind_rows(antigenic_map, antigenic_map_tmp2)
    }
}

## Randomize model parameters
par_tab <- generate_start_tab(par_tab)

## Simulate realistic-ish attack rates
attack_rates <- simulate_attack_rates(strain_isolation_times, 0.15,1,FALSE)
attack_rates[attack_rates>1] <- 1
plot(attack_rates)

## Simulate the data
sim_data <- simulate_data(par_tab=par_tab, group=1, n_indiv=n_indivs, 
                     strain_isolation_times=strain_isolation_times,
                     measured_strains=sampled_viruses,
                     sampling_times=sampling_times, nsamps=n_samps, 
                     antigenic_map=antigenic_map, 
                     titre_sensoring=0.2, ## Randomly censor 20% of measurements
                     age_min=age_min,age_max=age_max,
                     attack_rates=attack_rates, repeats=repeats,
                     mu_indices = NULL, 
                     measurement_indices = measurement_indices,
                     obs_dist=obs_type_dists)

sum(sim_data$infection_histories)
plot_data(sim_data$data,sim_data$infection_histories,strain_isolation_times,n_indivs=3) 

titre_dat <- sim_data$data %>% tidyr::drop_na()
titre_dat <- titre_dat %>% left_join(sim_data$ages)
titre_dat <- titre_dat %>% arrange(individual, obs_type,samples, virus, run) %>% as.data.frame()

## Save titre data
write_csv(titre_dat, file=paste0(save_wd,"/",run_name,"_titre_data.csv"))
## Save parameter table
write_csv(par_tab, file=paste0(save_wd,"/",run_name,"_par_tab.csv"))
## Save attack rates
write_csv(sim_data$attack_rates, file=paste0(save_wd,"/",run_name,"_attack_rates.csv"))
## Save infection histories
write_csv(as.data.frame(sim_data$infection_histories), file=paste0(save_wd,"/",run_name,"_infection_histories.csv"))


## Going to fit model to the data 3 times. Once with each observation type, and once with both.
data_type_list <- list(1, c(1,2))
run_names <- c("obs1","both")
obs_type_list <- list(1,c(1,2))

## Fit the model either with just obs_type 1, or both obs types
for(run in seq_along(data_type_list)){
#for(run in 2){
    run_name_use <- paste0(run_name, "_", run_names[run])
    obs_type_use <- obs_type_list[[run]]
    data_type_vector <- data_type_list[[run]]
    
    chain_wd_use <- paste0(chain_wd, "_",run_names[run])
    save_wd_use <- paste0(save_wd, "_",run_names[run])
    
    if(!dir.exists(save_wd_use)) dir.create(save_wd_use,recursive = TRUE)
    if(!dir.exists(chain_wd_use)) dir.create(chain_wd_use,recursive = TRUE)
    
    ## Only use data and parameter table relevant to these obs_types
    par_tab_use <- par_tab[par_tab$obs_type %in% obs_type_use,]
    antigenic_map_use <- antigenic_map[antigenic_map$obs_type %in% obs_type_use,]
    titre_dat_use <- titre_dat[titre_dat$obs_type %in% obs_type_use,]
    
    f <- create_posterior_func(par_tab_use,titre_dat_use,antigenic_map=antigenic_map_use,
                               strain_isolation_times = strain_isolation_times,
                               version=prior_version,solve_likelihood=TRUE,n_alive=NULL,
                               measurement_indices_by_time=measurement_indices,
                               data_type=data_type_vector
                               )

    ## Time runs and use dopar to run multiple chains in parallel
    t1 <- Sys.time()
    filenames <- paste0(chain_wd_use, "/",run_name_use, "_",1:n_chains)
    
    if(rerun){
    res <- foreach(x = filenames, .packages = c('data.table','plyr',"dplyr","tidyverse")) %dopar% {
        devtools::load_all(serosolver_wd)
        index <- 1
        lik <- -Inf
        inf_hist_correct <- 1
        while((!is.finite(lik) || inf_hist_correct > 0) & index < 100){
            start_tab <- generate_start_tab(par_tab_use)
            start_inf <- setup_infection_histories_total(titre_dat_use,strain_isolation_times,2,3)
            
            inf_hist_correct <- sum(check_inf_hist(titre_dat_use, strain_isolation_times, start_inf))
            y <- f(start_tab$values, start_inf)
            lik <- sum(y[[1]])
            index <- index + 1
        }
        
        write.csv(start_tab, paste0(x, "_start_tab.csv"))
        write.csv(start_inf, paste0(x, "_start_inf_hist.csv"))
        write.csv(antigenic_map, paste0(x, "_antigenic_map.csv"))
        write.csv(titre_dat, paste0(x, "_titre_dat.csv"))
        
        res <- serosolver::run_MCMC(start_tab, 
                                    titre_dat_use, 
                                    antigenic_map_use, 
                                    start_inf_hist=start_inf,
                                    filename=x,
                                    CREATE_POSTERIOR_FUNC=create_posterior_func, 
                                    CREATE_PRIOR_FUNC = NULL,
                                    version=prior_version,
                                    mcmc_pars=mcmc_pars,
                                    measurement_indices= measurement_indices, ## measurement_indices, ## NULL
                                    measurement_random_effects = FALSE, ## TRUE, ## FALSE
                                    solve_likelihood=TRUE,
                                    data_type=data_type_vector)
    }
    }
    run_time_fast <- Sys.time() - t1
    run_time_fast
    
    ## Read in chains for trace plot
    chains <- load_mcmc_chains(chain_wd_use,convert_mcmc=TRUE,burnin = mcmc_pars["adaptive_period"],unfixed = TRUE)
    pdf(paste0(save_wd_use,"/",run_name_use,"_chain.pdf"))
    plot(as.mcmc.list(chains$theta_list_chains))
    dev.off()
    
    ## Plot posterior densities and compare to ground truth parameters
    theta_chain <- chains$theta_chain
    theta_chain_melted <- reshape2::melt(as.data.frame(theta_chain),id.vars=c("sampno","chain_no")) %>% as_tibble()
    
    par_tab_tmp <-  par_tab_use %>% group_by(names) %>% mutate(n=1:n()) %>% mutate(variable = ifelse(n > 1, paste0(names,".",n-1), names)) %>%
        select(-n) %>% mutate(obs_type_p = paste0("Obs type: ", obs_type))
    theta_chain_melted <- theta_chain_melted %>% left_join(par_tab_tmp)%>% mutate(obs_type_p = paste0("Obs type: ", obs_type))
    
   
    p_estimates_theta <- ggplot(theta_chain_melted %>% filter(fixed == 0, type == 1)) +
        geom_density(aes(x=value,fill=obs_type_p),alpha=0.4) + 
        geom_vline(data=par_tab_tmp%>% filter(fixed == 0, type == 1), 
                   aes(xintercept=values,col="True value"),linetype="dashed") +
        facet_wrap(names~obs_type_p,scales="free") +
        scale_color_manual(name="Estimate: ",values=c("True value"="grey40"))+
        scale_fill_manual(name="Observation type: ",values=c("Obs type: 1"="blue","Obs type: 2"="red")) +
        theme_minimal() +
        ylab("Density") +
        xlab("Value") +
        theme(legend.position="bottom")

    p_estimates_bias <- ggplot(theta_chain_melted %>% filter(fixed == 0, type == 3)) +
        geom_density(aes(x=value,fill=obs_type_p),alpha=0.4) + 
        geom_vline(data=par_tab_tmp%>% filter(fixed == 0, type == 3), 
                   aes(xintercept=values,col="True value"),linetype="dashed") +
        facet_wrap(variable~obs_type_p,scales="free") +
        geom_vline(xintercept=0) +
        scale_x_continuous(limits=c(-1.25,1.25),breaks=seq(-1.5,1.5,by=0.5)) +
        scale_color_manual(name="Estimate: ",values=c("True value"="grey40"))+
        scale_fill_manual(name="Observation type: ",values=c("Obs type: 1"="blue","Obs type: 2"="red")) +
        theme_minimal() +
        ylab("Density") +
        xlab("Value") +
        theme(legend.position="bottom")

    ggsave(paste0(save_wd_use,"/",run_name_use,"_par_estimates.pdf"),p_estimates_theta,height=8,width=7,units="in",dpi=300)
    ggsave(paste0(save_wd_use,"/",run_name_use,"_rho_estimates.pdf"),p_estimates_bias,height=6,width=7,units="in",dpi=300)
    
    
    ## Read in chains for all other plots
    chains <- load_mcmc_chains(chain_wd_use,convert_mcmc=FALSE,burnin = mcmc_pars["adaptive_period"],unfixed = FALSE)
    chain <- as.data.frame(chains$theta_chain)
    inf_chain <- chains$inf_chain 
    
    ## Plot kinetics parameter estimates and number of infections
    p1 <- plot_posteriors_theta(chain, par_tab, burnin=0,samples=50,TRUE,FALSE,FALSE, TRUE,"")
    p_total_inf <- plot_total_number_infections(inf_chain,FALSE)
    ggsave(paste0(save_wd_use,"/",run_name_use,"_total_inf.pdf"),p_total_inf[[1]],height=5,width=7,units="in",dpi=300)
    
    p_cumu_infs <- generate_cumulative_inf_plots(inf_chain,indivs=1:5,real_inf_hist=sim_data$infection_histories,strain_isolation_times = strain_isolation_times,nsamp=50)
    ggsave(paste0(save_wd_use,"/",run_name_use,"_inf_hists.pdf"),p_cumu_infs[[1]],height=8,width=7,units="in",dpi=300)
    
    n_alive <- get_n_alive_group(titre_dat,strain_isolation_times)
    
    ## Plot attack rates
    n_inf <- sim_data$infection_histories %>% colSums()
    true_ar <- data.frame(j=strain_isolation_times,group=1,AR=n_inf/n_alive[1,])
    p_ar <- plot_attack_rates(inf_chain, titre_dat_use, 2000:2013,
                              pad_chain=FALSE,plot_den=TRUE,n_alive = n_alive,
                              true_ar=true_ar,
                              prior_pars = c("prior_version"=2,"alpha"=1/3,"beta"=1/3),
                              by_val=1)
    p_ar
    ggsave(paste0(save_wd_use,"/",run_name_use,"_ar.pdf"),p_ar,height=7,width=8,units="in",dpi=300)
    
    
    ## Plot model fits, different subsets depending on which observation types we're using
    titre_preds <- get_titre_predictions(chain = chain[chain$chain_no == 1,], 
                                         infection_histories = inf_chain[inf_chain$chain_no == 1,], 
                                         titre_dat = titre_dat_use, 
                                         individuals = unique(titre_dat_use$individual),
                                         antigenic_map = antigenic_map_use, 
                                         par_tab = par_tab_use,expand_titredat=FALSE,
                                         nsamp=500,
                                         data_type = data_type_vector)
    
    use_indivs <- sample(unique(titre_dat_use$individual), 9)
    
    
    titre_pred_p <- ggplot(titre_preds$predictions[titre_preds$predictions$individual %in% use_indivs,])+
        geom_ribbon(data=titre_preds$predicted_observations[titre_preds$predicted_observations$individual %in% use_indivs,], 
                    aes(x=virus,ymin=lower,ymax=upper,fill=as.factor(obs_type)),alpha=0.1) +
        geom_ribbon(aes(x=virus,ymin=lower, ymax=upper,fill=as.factor(obs_type)),alpha=0.25)+
        geom_line(aes(x=virus, y=median,col=as.factor(obs_type)))+
        geom_point(aes(x=virus, y=titre,col=as.factor(obs_type)))+
        coord_cartesian(ylim=c(0,8))+
        ylab("log titre") +
        xlab("Time of virus circulation") +
        theme_classic() +
        facet_wrap(obs_type~individual)
    
    ## Model fits for observation type 1
    titre_pred_p1 <- ggplot(titre_preds$predictions[titre_preds$predictions$individual %in% use_indivs,] %>% filter(obs_type == 1))+
        geom_ribbon(data=titre_preds$predicted_observations[titre_preds$predicted_observations$individual %in% use_indivs,]%>% filter(obs_type == 1), 
                    aes(x=virus,ymin=lower,ymax=upper),alpha=0.1,fill="blue") +
        geom_ribbon(aes(x=virus,ymin=lower, ymax=upper),alpha=0.25,fill="blue")+
        geom_line(aes(x=virus, y=median),col="blue")+
        geom_point(data=titre_dat_use[titre_dat_use$individual %in% use_indivs,] %>% filter(obs_type==1), aes(x=virus, y=titre),col="blue")+
        coord_cartesian(ylim=c(0,8))+
        ylab("log titre") +
        xlab("Time of virus circulation") +
        ggtitle("Model fits to observation data type 1") +
        theme_classic() +
        facet_grid(individual~samples)
    
    titre_pred_p1
    ggsave(paste0(save_wd_use,"/",run_name_use,"_obs_type_1_fits.pdf"),titre_pred_p1,height=7,width=8,units="in",dpi=300)
    
    ## Model fits for observation type 2
    if(2 %in% obs_type_use){
    titre_pred_p2 <- ggplot(titre_preds$predictions[titre_preds$predictions$individual %in% use_indivs,] %>% filter(obs_type == 2))+
        geom_ribbon(data=titre_preds$predicted_observations[titre_preds$predicted_observations$individual %in% use_indivs,]%>% filter(obs_type == 2), 
                    aes(x=virus,ymin=lower,ymax=upper),alpha=0.1,fill="red") +
        geom_ribbon(aes(x=virus,ymin=lower, ymax=upper),alpha=0.25,fill="red")+
        geom_line(aes(x=virus, y=median),col="red")+
        geom_point(data=titre_dat_use[titre_dat_use$individual %in% use_indivs,] %>% filter(obs_type==2),aes(x=virus, y=titre),col="red")+
        coord_cartesian(ylim=c(0,8))+
        ylab("log titre") +
        xlab("Time of virus circulation") +
        ggtitle("Model fits to observation data type 2") +
        theme_classic() +
        facet_grid(individual~samples)
    
        titre_pred_p2
        ggsave(paste0(save_wd_use,"/",run_name_use,"_obs_type_2_fits.pdf"),titre_pred_p2,height=7,width=8,units="in",dpi=300)
    }
    if((1 %in% obs_type_use) & (2 %in% obs_type_use)){
        titre_pred_combined <- ggplot(titre_preds$predictions[titre_preds$predictions$individual %in% use_indivs,] )+
            #geom_ribbon(data=titre_preds$predicted_observations[titre_preds$predicted_observations$individual %in% use_indivs,], 
            #            aes(x=virus,ymin=lower,ymax=upper,fill=as.factor(obs_type)),alpha=0.1) +
            geom_ribbon(aes(x=virus,ymin=lower, ymax=upper,fill=as.factor(obs_type)),alpha=0.25)+
            geom_line(aes(x=virus, y=median,col=as.factor(obs_type)))+
            geom_point(data=titre_dat_use[titre_dat_use$individual %in% use_indivs,], aes(x=virus, y=titre,col=as.factor(obs_type)),alpha=0.5)+
            coord_cartesian(ylim=c(0,8))+
            ylab("log titre") +
            xlab("Time of virus circulation") +
            scale_color_manual(name="Observation type",values=c("1"="blue","2"="red")) +
            scale_fill_manual(name="Observation type",values=c("1"="blue","2"="red")) +
            theme_classic() +
            facet_grid(individual~samples)
        ggsave(paste0(save_wd_use,"/",run_name_use,"_both_fits.pdf"),titre_pred_combined,height=7,width=8,units="in",dpi=300)
        
    }
}
