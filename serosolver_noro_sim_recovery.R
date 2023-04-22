######################################################
## SIMULATION RECOVERY TEST -- NOROVIRUS DATA
## Author: James Hay
## Date: 03 March 2023
## Summary: simulates some serosurvey data and fits serosolver

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
#library(devtools)

library(misc3d)
library(plot3D)
library(rgl)

#serosolver_wd <- "~/Documents/GitHub/serosolver/"
#devtools::load_all(serosolver_wd)
library(serosolver)
packageDescription("serosolver") # check it is the right branch

# notes: a "makefile" is needed to run this. I think this si something to do with ensuring the right R packages 
# are installed.

rm(list = ls())

run_name <- "data_test_ac" #"sim_noro_allsampled" #
main_wd <-  "/Users/lsh1603970/GitHub/serosolver_norovirus_2sera" #"~/Documents/norovirus_test/"
setwd(main_wd)
chain_wd <- paste0(main_wd,"/chains/",run_name)
save_wd <- paste0(main_wd,"/figures/chain_plots/")

if(!dir.exists(save_wd)) dir.create(save_wd,recursive = TRUE)
if(!dir.exists(chain_wd)) dir.create(chain_wd,recursive = TRUE)

buckets <- 1 ## Ignore
prior_version <- 2 ## Which prior on the infection histories to use? Prior version 2 is generally preferred
n_chains <- 3 ## Number of MCMC chains to run

rerun <- TRUE ## Set to FALSE if you just want to load in previously run chains

setwd(main_wd)
print(paste0("In directory: ", main_wd))
print(paste0("Saving to: ", save_wd))

set.seed(1)

## Run multiple chains in parallel
cl <- makeCluster(n_chains)
registerDoParallel(cl)

## MCMC settings, not super important but can be tweaked
mcmc_pars <- c("save_block"=100,"thin"=100,"thin_hist"=500,"iterations"=100000,
               "adaptive_period"=50000,
               "burnin"=0,"switch_sample"=2,"hist_switch_prob"=0.05,
               "year_swap_propn"=0.8,"swap_propn"=0.5,
               "inf_propn"=0.5,"hist_sample_prob"=1,"move_size"=3, "hist_opt"=0,
               "popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,propose_from_prior=TRUE)

## Simulation parameters
n_indivs <- 250
n_groups <- 1 ## Leave as 1, can in theory make multiple groups with distinct FOI parameters, but needs extra setup
n_samps <- 1 ## Number of samples per person
repeats <- 1 ## Number of repeat measurements per variant/sample combination
samp_min <- 2008 ## First sample year
samp_max <- 2012 ## Final sample year
year_min <- 2000 ## First year of possible circulation (ie. time 0 of the simulation)
year_max <- 2012 ## Final year of possible circulation
age_min <- 1 ## Age minimum and maximum in years, simulated from a uniform distribution
age_max <- 7

## Viruses and times for samples
sampled_viruses <- c(2002,2006,2009,2012)
sampling_times <- seq(samp_min, samp_max, by=1)

if(run_name == "data_test_ac"){
  antigenic_mapB <- read_csv("antigenic_map_noro_inferred_jump.csv")
  oo <- match(sampled_viruses,antigenic_mapB$inf_times)
  antigenic_coords <- antigenic_mapB[oo,]
  names(antigenic_coords) <- c("X","Y","Strain")
  antigenic_map <- (read_csv("antigenic_map_noro_inferred_trueequal.csv")) #
}
if(run_name == "data_test_ac2"){
  antigenic_mapB <- read_csv("antigenic_map_noro_inferred_jump.csv")
  oo <- match(sampled_viruses,antigenic_mapB$inf_times)
  antigenic_coords <- antigenic_mapB[oo,]
  names(antigenic_coords) <- c("X","Y","Strain")
  antigenic_map <- as.data.frame(read_csv("antigenic_map_noro_inferred_temporal.csv")) #
}
if(run_name != "data_test_ac2" & run_name != "data_test_ac"){
  ## Create a fake antigenic map -- can put in what you like here
  antigenic_coords <- data.frame(Strain=c(2000,2002,2006,2009,2012),X=c(0,0.5,3,3.5,4),Y=c(0,2,1,3,4))
  antigenic_map <- generate_antigenic_map_flexible(antigenic_coords,
                                                   year_max=2013,year_min=2000,spar = 0.001)
}

ggplot(data=antigenic_map,aes(x=x_coord,y=y_coord)) + geom_line() +
  geom_text(aes(x=x_coord+0.5,y=y_coord,label=inf_times)) +
  geom_point(data=antigenic_coords,aes(x=X,y=Y),col="blue") +
  geom_text(data=antigenic_coords,aes(x=X+0.5,y=Y,label=Strain),col="blue")

ggplot(antigenic_map,aes(x=inf_times,y=x_coord)) + geom_point()

#ggplot()

strain_isolation_times <- antigenic_map$inf_times
n_times <- length(strain_isolation_times)

## Set up parameter table
par_tab <- read.csv("par_tab_base_new.csv",stringsAsFactors=FALSE)
par_tab <- par_tab[par_tab$names != "phi",]
par_tab[par_tab$names == c("alpha","beta"),c("values")] <- c(1,2) ## Can also try c(1,1), or something informative. Just the parameters of a beta distribution which acts as the prior on the per-time attack rate.

#hist(rbeta(1000,1,2))

## Just some setup of the parameter table to get parameters vaguely like you showed me
par_tab$fixed <- 1
par_tab[par_tab$names %in% c("mu","mu_short","wane","sigma1","sigma2","tau","obs_sd"),"fixed"] <- 0
par_tab[par_tab$names %in% c("mu","mu_short","wane","sigma1","sigma2","tau","obs_sd"),"values"] <- c(2,5,0.4,0.1,0.01,0.5,1)


## Create some random measurement offsets -- these add a systematic shift to each observed variant. Only those sampled variant years have offsets, the rest are 0
## These are unknown parameters to be estimated, assuming the offsets are drawn from ~ norm(0, 1)
par_tab_rhos <- data.frame(names="rho",values=rep(0,n_times),fixed=1,steps=0.1,
                           lower_bound=-3,upper_bound=3,lower_start=-1,
                           upper_start=1,type=3)
par_tab_rhos[which(strain_isolation_times%in%sampled_viruses),"values"] <- rnorm(length(sampled_viruses), 1)
par_tab_rhos[which(strain_isolation_times%in%sampled_viruses),"fixed"] <-0
measurement_indices <- seq_along(strain_isolation_times)
par_tab <- bind_rows(par_tab, par_tab_rhos)

## Simulate realistic-ish attack rates
attack_rates <- simulate_attack_rates(strain_isolation_times, 0.15,1,FALSE)
attack_rates[attack_rates>1] <- 1
plot(strain_isolation_times,attack_rates,pch=19)

## Simulate the data
sim_data <- simulate_data(par_tab=par_tab, group=1, n_indiv=n_indivs, 
                     #buckets=buckets,
                     strain_isolation_times=strain_isolation_times,
                     measured_strains=sampled_viruses,
                     sampling_times=sampling_times, nsamps=n_samps, 
                     antigenic_map=antigenic_map, 
                     titre_sensoring=0.2, ## Randomly censor 20% of measurements
                     age_min=age_min,age_max=age_max,
                     attack_rates=attack_rates, repeats=repeats,
                     mu_indices = NULL, measurement_indices = measurement_indices)#,
                    # add_noise=TRUE
                    #)

sum(sim_data$infection_histories)
plot_data(sim_data$data,sim_data$infection_histories,strain_isolation_times,n_indivs=10)

titre_dat <- sim_data$data %>% tidyr::drop_na()
titre_dat <- titre_dat %>% group_by(individual, samples, virus) %>% mutate(run=1:n()) %>% ungroup()
titre_dat <- titre_dat %>% left_join(sim_data$ages)
titre_dat <- titre_dat %>% arrange(individual, samples, virus, run) %>% as.data.frame()

ggplot(titre_dat,aes(x=samples,y=titre,group=samples)) + geom_boxplot() + facet_wrap(~virus)


## Save titre data
write_csv(titre_dat, file=paste0(save_wd,"/",run_name,"_titre_data.csv"))
## Save parameter table
write_csv(par_tab, file=paste0(save_wd,"/",run_name,"_par_tab.csv"))
## Save attack rates
write_csv(sim_data$attack_rates, file=paste0(save_wd,"/",run_name,"_attack_rates.csv"))
## Save infection histories
write_csv(as.data.frame(sim_data$infection_histories), file=paste0(save_wd,"/",run_name,"_infection_histories.csv"))

head(titre_dat)
table(samp=titre_dat$samples,virus=titre_dat$virus)

f <- create_posterior_func(par_tab,titre_dat,antigenic_map=antigenic_map,strain_isolation_times = strain_isolation_times,
                           version=prior_version,solve_likelihood=TRUE,n_alive=NULL,
                           measurement_indices_by_time=measurement_indices ## NULL
                           )

## Time runs and use dopar to run multiple chains in parallel
#library(data.table)

# check which parameters we are estimating
par_tab[par_tab$fixed==0,]

t1 <- Sys.time()
filenames <- paste0(chain_wd, "/",run_name, "_",1:n_chains)
if(rerun){
  res <- foreach(x = filenames, .packages = c('serosolver','data.table','plyr',"dplyr")) %dopar% {
      devtools::load_all(save_wd)
      index <- 1
      lik <- -Inf
      inf_hist_correct <- 1
      while((!is.finite(lik) || inf_hist_correct > 0) & index < 100){
          start_tab <- generate_start_tab(par_tab)
          start_inf <- setup_infection_histories_total(titre_dat,strain_isolation_times,2,3)
          
          inf_hist_correct <- sum(check_inf_hist(titre_dat, strain_isolation_times, start_inf))
          y <- f(start_tab$values, start_inf)
          lik <- sum(y[[1]])
          index <- index + 1
      }
      write.csv(start_tab, paste0(x, "_start_tab.csv"))
      write.csv(start_inf, paste0(x, "_start_inf_hist.csv"))
      write.csv(antigenic_map, paste0(x, "_antigenic_map.csv"))
      write.csv(titre_dat, paste0(x, "_titre_dat.csv"))
      res <- serosolver::run_MCMC(start_tab, titre_dat, antigenic_map, 
                                  start_inf_hist=start_inf,filename=x,
                                  CREATE_POSTERIOR_FUNC=create_posterior_func, 
                                  CREATE_PRIOR_FUNC = NULL,
                                  version=prior_version,
                                  mcmc_pars=mcmc_pars,
                                  measurement_indices=measurement_indices, ## NULL
                                  measurement_random_effects = TRUE, ## FALSE
                                  solve_likelihood=TRUE)
  }
}
run_time_fast <- Sys.time() - t1
run_time_fast # 20 mins. 

## Read in chains for trace plot
chains <- load_mcmc_chains(chain_wd,convert_mcmc=TRUE,burnin = mcmc_pars["adaptive_period"],unfixed = TRUE)
pdf(paste0(save_wd,"/",run_name,"_chain.pdf"))
plot(as.mcmc.list(chains$theta_list_chains))
dev.off()

## Read in chains for all other plots
chains <- load_mcmc_chains(chain_wd,convert_mcmc=FALSE,burnin = mcmc_pars["adaptive_period"],unfixed = FALSE)
chain <- as.data.frame(chains$theta_chain)
inf_chain <- chains$inf_chain 

#summary stats
library(coda)
list_chains1 <- chains[[1]]
tmp <- summary(chain[,c("mu","sigma1",
                         "tau",
                         #"wane",
                         "error","total_infections",
                         "lnlike","prior_prob")])
tmp

# how well do we re-estimate parameters?
df2 <- data.frame(variable=c("mu","mu_short","tau","sigma1","sigma2","wane","error"),
                  value=c(2,5,0.4,0.01,0.5,0.1,1))
p1 <- list_chains1 %>% select(mu,mu_short,tau,sigma1,sigma2,wane,error) %>% tidyr::gather(variable, value) %>%
  ggplot(aes(value)) + geom_histogram(bins = 51) +
  geom_vline(data=df2,aes(xintercept=value),col="red") +
  facet_wrap(~variable, scales = 'free_x') 

p1

## Plot kinetics parameter estimates and number of infections
p1 <- plot_posteriors_theta(chain, par_tab, burnin=0,samples=100,TRUE,FALSE,FALSE, TRUE,"")
p_total_inf <- plot_total_number_infections(inf_chain,FALSE)
ggsave(paste0(save_wd,"/",run_name,"_total_inf.pdf"),p_total_inf[[1]],height=5,width=7,units="in",dpi=300)

p_cumu_infs <- generate_cumulative_inf_plots(inf_chain,indivs=1:5,real_inf_hist=sim_data$infection_histories,strain_isolation_times = strain_isolation_times)
ggsave(paste0(save_wd,"/",run_name,"_inf_hists.pdf"),p_cumu_infs[[1]],height=8,width=7,units="in",dpi=300)

n_alive <- get_n_alive_group(titre_dat,strain_isolation_times)

## Plot attack rates
n_inf <- sim_data$infection_histories %>% colSums()
true_ar <- data.frame(j=strain_isolation_times,group=1,AR=n_inf/n_alive[1,])
p_ar <- plot_attack_rates(inf_chain, titre_dat, 2000:2013,
                          pad_chain=FALSE,plot_den=TRUE,n_alive = n_alive,
                          true_ar=true_ar,
                          prior_pars = c("prior_version"=2,"alpha"=1,"beta"=1))
p_ar
ggsave(paste0(save_wd,"/",run_name,"_ar.pdf"),p_ar,height=7,width=8,units="in",dpi=300)


## Plot model fits
titre_preds <- get_titre_predictions(chain = chain[chain$chain_no == 1,], 
                                     infection_histories = inf_chain[inf_chain$chain_no == 1,], 
                                     titre_dat = titre_dat, 
                                     individuals = unique(titre_dat$individual),
                                     antigenic_map = antigenic_map, 
                                     par_tab = par_tab,expand_titredat=FALSE)

use_indivs <- sample(unique(titre_dat$individual), 9)


titre_pred_p <- ggplot(titre_preds$predictions[titre_preds$predictions$individual %in% use_indivs,])+
    geom_ribbon(data=titre_preds$predicted_observations[titre_preds$predicted_observations$individual %in% use_indivs,], 
                aes(x=virus,ymin=lower,ymax=upper),fill="grey90") +
    geom_ribbon(aes(x=virus,ymin=lower, ymax=upper),fill="gray70")+
    geom_line(aes(x=virus, y=median))+
    geom_point(aes(x=virus, y=titre))+
    coord_cartesian(ylim=c(0,8))+
    ylab("log titre") +
    xlab("Time of virus circulation") +
    theme_classic() +
    facet_wrap(~individual)
titre_pred_p
ggsave(paste0(save_wd,"/",run_name,"_titre_fits.pdf"),titre_pred_p,height=7,width=8,units="in",dpi=300)

