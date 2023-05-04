# plot model outputs

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

run_name <- "data_t2"
main_wd <- "~/GitHub/serosolver_norovirus_2sera/"
setwd(main_wd)

chain_wd <- paste0(main_wd,"/chains_avidity/",run_name)
save_wd <- paste0(main_wd,"/figures/chain_plots_avidity/") # where the data are?

# assumign models have been run in "ss_noro_data_avidity.r"

titre_dat <- read_csv(file=paste0(save_wd,"/",run_name,"_titre_data.csv"))

# MCMC settings, not super important but can be tweaked
mcmc_pars <- c("save_block"=100,"thin"=10,"thin_hist"=50,
               "iterations"=200000,
               "adaptive_period"=50000,
               "burnin"=0,"switch_sample"=2,"hist_switch_prob"=0.05,
               "year_swap_propn"=0.8,"swap_propn"=0.5,
               "inf_propn"=0.5,"hist_sample_prob"=1,"move_size"=3, "hist_opt"=0,
               "popt"=0.44,"popt_hist"=0.44,"opt_freq"=2000,propose_from_prior=TRUE)
sampled_viruses <- c(2002,2006,2009,2012)
n_obs_types <- 2
## Set up parameter table
par_tab <- read.csv("par_tab_base_new.csv",stringsAsFactors=FALSE)
par_tab <- par_tab[par_tab$names != "phi",]
par_tab[par_tab$names %in% c("alpha","beta"),c("values")] <- c(1,2) ## Can also try c(1,1), or something informative. Just the parameters of a beta distribution which acts as the prior on the per-time attack rate.

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

# plots
run_names <- c("obs1","both")
run <- 2
run_name_use <- paste0(run_name, "_", run_names[run])
par_tab_use <- par_tab[par_tab$obs_type %in% obs_type_use,]

chain_wd_use <- paste0(chain_wd, "_",run_names[run])
save_wd_use <- paste0(save_wd, "_",run_name[run])

## Read in chains for trace plot
chains <- load_mcmc_chains(chain_wd_use,convert_mcmc=TRUE,burnin = mcmc_pars["adaptive_period"],unfixed = TRUE)
pdf(paste0(save_wd_use,"/",run_name_use,"_chain.pdf"))
plot(as.mcmc.list(chains$theta_list_chains))
dev.off()

## Plot posterior densities and compare to ground truth parameters
theta_chain <- chains$theta_chain
theta_chain_melted <- reshape2::melt(as.data.frame(theta_chain),id.vars=c("sampno","chain_no")) %>% as_tibble()

head(theta_chain_melted)
head(theta_chain)

# can I plot correlation here?
library(cowplot)
tmp <- as.data.frame(theta_chain) %>% filter(chain_no == 1)

# mega correlation estimates
vals <- names(tmp)
vals_use <- vals[c(c(5:10),c(21:26))]
ss <- sample(dim(tmp)[1],1000,replace = FALSE)
cor_dat <- tmp[ss,vals_use]

p1 <- ggplot(cor_dat,aes(x=mu_short,y=mu_short.1)) + geom_point()
p2 <- ggplot(cor_dat,aes(x=mu_short,y=tau.1)) + geom_point()
p3 <- ggplot(cor_dat,aes(x=mu_short,y=tau.1)) + geom_point()
p4 <- ggplot(cor_dat,aes(x=sigma1,y=sigma2)) + geom_point()

pdf(paste0(run_name_use,"_correlation_4_plot.pdf"),height=3,width=9)
plot_grid(p1,p2,p3,p4,ncol=4)
dev.off()

library(Hmisc)
res <- rcorr(as.matrix(cor_dat))
round(res$r, 2)
round(res$P, 5)
library("PerformanceAnalytics")
chart.Correlation(cor_dat, histogram=TRUE, pch=20)

pdf(paste0(run_name_use,"_correlation_plot.pdf"),height=8,width=9)
chart.Correlation(cor_dat, histogram=TRUE, pch=20)
dev.off()

#

par_tab_tmp <-  par_tab_use %>% group_by(names) %>% mutate(n=1:n()) %>% mutate(variable = ifelse(n > 1, paste0(names,".",n-1), names)) %>%
  select(-n) %>% mutate(obs_type_p = paste0("Obs type: ", obs_type))
theta_chain_melted <- theta_chain_melted %>% left_join(par_tab_tmp)%>% mutate(obs_type_p = paste0("Obs type: ", obs_type))

p_estimates_theta <- ggplot(theta_chain_melted %>% filter(fixed == 0, type == 1)) +
  geom_density(aes(x=value,fill=obs_type_p),alpha=0.4) + 
  #geom_vline(data=par_tab_tmp%>% filter(fixed == 0, type == 1), 
  #           aes(xintercept=values,col="True value"),linetype="dashed") +
  facet_wrap(names~obs_type_p,scales="free") +
  #scale_color_manual(name="Estimate: ",values=c("True value"="grey40"))+
  scale_fill_manual(name="Observation type: ",values=c("Obs type: 1"="blue","Obs type: 2"="red")) +
  theme_minimal() +
  ylab("Density") +
  xlab("Value") +
  theme(legend.position="bottom")

ggsave(paste0(save_wd_use,"/",run_name_use,"_par_estimates.pdf"),p_estimates_theta,height=8,width=7,units="in",dpi=300)

p_estimates_bias <- ggplot(theta_chain_melted %>% filter(fixed == 0, type == 3)) +
  geom_density(aes(x=value,fill=obs_type_p),alpha=0.4) + 
  #geom_vline(data=par_tab_tmp%>% filter(fixed == 0, type == 3), 
  #           aes(xintercept=values,col="True value"),linetype="dashed") +
  facet_wrap(variable~obs_type_p,scales="free") +
  geom_vline(xintercept=0) +
  scale_x_continuous(limits=c(-1.25,1.25),breaks=seq(-1.5,1.5,by=0.5)) +
  #scale_color_manual(name="Estimate: ",values=c("True value"="grey40"))+
  scale_fill_manual(name="Observation type: ",values=c("Obs type: 1"="blue","Obs type: 2"="red")) +
  theme_minimal() +
  ylab("Density") +
  xlab("Value") +
  theme(legend.position="bottom")

ggsave(paste0(save_wd_use,"/",run_name_use,"_rho_estimates.pdf"),p_estimates_bias,height=6,width=7,units="in",dpi=300)

## Read in chains for all other plots
chains <- load_mcmc_chains(chain_wd_use,convert_mcmc=FALSE,burnin = mcmc_pars["adaptive_period"],unfixed = FALSE)
chain <- as.data.frame(chains$theta_chain)
inf_chain <- chains$inf_chain 

## Plot kinetics parameter estimates and number of infections
p1 <- plot_posteriors_theta(chain, par_tab, burnin=0,samples=50,TRUE,FALSE,FALSE, TRUE,"")
p_total_inf <- plot_total_number_infections(inf_chain,FALSE)
ggsave(paste0(save_wd_use,"/",run_name_use,"_total_inf.pdf"),p_total_inf[[1]],height=5,width=7,units="in",dpi=300)

p_cumu_infs <- generate_cumulative_inf_plots(inf_chain,indivs=1:5,
                                             #real_inf_hist=sim_data$infection_histories,
                                             strain_isolation_times = strain_isolation_times,nsamp=50)
ggsave(paste0(save_wd_use,"/",run_name_use,"_inf_hists.pdf"),p_cumu_infs[[1]],height=8,width=7,units="in",dpi=300)

n_alive <- get_n_alive_group(titre_dat,strain_isolation_times)

## Plot attack rates
#n_inf <- sim_data$infection_histories %>% colSums()
#true_ar <- data.frame(j=strain_isolation_times,group=1,AR=n_inf/n_alive[1,])
p_ar <- plot_attack_rates(inf_chain, titre_dat_use, 2000:2013,
                          pad_chain=FALSE,plot_den=TRUE,n_alive = n_alive,
                          #true_ar=true_ar,
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
    facet_wrap(individual)
  
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
    facet_wrap(~individual)
  ggsave(paste0(save_wd_use,"/",run_name_use,"_both_fits.pdf"),titre_pred_combined,height=7,width=8,units="in",dpi=300)
  
}

# end