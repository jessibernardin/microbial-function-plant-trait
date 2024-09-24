#### Microbial function mediates leaf traits in a pitcher plant model system ####
#### Authors: Jessica R. Bernardin, Erica B. Young, Leonora S. Bittleston ####
#### last update : July 30, 2024 ####
#### Physiological Functions Analysis

#### LOAD REQUIRED PACKAGES ####
packages_to_load <- c(
  "ggplot2", "vegan", "lme4", "tidyverse", "effects", "growthcurver",
  "plyr", "dplyr", "reshape", "reshape2", "ape", "DiagrammeR", "tidyverse",
  "tidybayes", "coefplot", "standardize", "bayesplot", "MCMCvis", "car",
  "patchwork", "ggpubr", "corrr", "ggcorrplot", "factoextra", "MASS",
  "pairwiseAdonis", "plotrix", "gridExtra", "multcompView", "ggeffects", "this.path", "brms", "ggmulti"
)

# Load and install required packages
for (i in packages_to_load) { #Installs packages if not yet installed
  if (!require(i, character.only = TRUE)) install.packages(i)
}

setwd(this.path::here())

#### READ IN THE DATA ####
#bacterial and plant metadata
data <- read.csv("micro_function_plant_trait_metadata.csv", header = TRUE, check.names = FALSE)

#filter out plants that had damage
data_filt <- data[!(data$plant_number %in% c("25", "44", "51", "9", "29", "NEG")), ]
data_filt$treatment <- as.factor(data_filt$treatment)
data_filt$day <- as.factor(data_filt$day)
data_filt$bacterial_respiration_rate_at_24hr <- as.numeric(data_filt$bacterial_respiration_rate_at_24hr)

#ecoplate data
ecodat <- read.csv("eco_all_2021.csv", header = TRUE)

#original cultures (day0) data
orig_resp <- read.csv("original_cultures_respiration.csv", header = TRUE)
orig_enz_rkn <- read.csv("original_cultures_enzyme_growth_data.csv", header = TRUE)

#colony forming units for inoculating cultures
m01v1 <- c(257, 255)
m06v1 <- c(194, 216)
m09v1 <- c(175, 216)

treatment <- rep(c("m01v1", "m06v1", "m09v1"), each = 2)
rep <- rep(c("rep1", "rep2"), times = 3)
value <- c(257, 255, 194, 216, 175, 216)

cfu <- data.frame(treatment = treatment, rep = rep, value = value)
aov_cfu <- aov(value ~ treatment, data = cfu)
summary(aov_cfu)
mean(cfu$value)

#### Microbial Functions of cultures prior to inoculation in planta ####
#growth rates of inoculating cultures
glm.growth_cult<- brm(growth_r ~ treatment, data= orig_enz_rkn, family=Gamma(link="log"), iter=10000)
summary(glm.growth_cult)
pp_check(glm.growth_cult)

#chitinase for inoculating cultures
glm.chit.org<- brm(chitinase_7_29 ~ treatment, data= orig_enz_rkn, family=Gamma(link="log"), iter=10000)
summary(glm.chit.org)
pp_check(glm.chit.org)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(orig_enz_rkn, varname="chitinase_7_29", 
                    groupnames=c("treatment"))
df2$treatment=as.factor(df2$treatment)

#Fig2A
ggplot(df2, aes(x = treatment, y = chitinase_7_29, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE, alpha=.75) +
  scale_fill_manual(values=c("#014b7a", "#ffaf49", "#44b7c2")) +theme_classic()+ ylab("Chitinase Activity")+
  geom_errorbar(aes(ymin=chitinase_7_29-sd, ymax=chitinase_7_29+sd), width=.2,
                position=position_dodge(.9))+theme(
    axis.text.x = element_text(color = "black", size = 16),
    axis.text.y = element_text(color = "black", size = 16),
    axis.title.x = element_text(color = "black", size = 16),
    axis.title.y = element_text(color = "black", size = 16)
  )

#protease of inoculating cultures
glm.prot.org<- brm(protease_7_29 ~ treatment, data= orig_enz_rkn, iter=10000)
#I did a normal distribution instead of gamma bc the pp_check was a better fit with normal
summary(glm.prot.org)
pp_check(glm.prot.org)
df2 <- data_summary(orig_enz_rkn, varname="protease_7_29", 
                    groupnames=c("treatment"))
df2$treatment=as.factor(df2$treatment)

#Fig2B
ggplot(df2, aes(x = treatment, y = protease_7_29, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE, alpha=0.75) +
  scale_fill_manual(values=c("#014b7a", "#ffaf49", "#44b7c2")) +theme_classic()+ ylab("Protease Activity")+
  geom_errorbar(aes(ymin=protease_7_29-sd, ymax=protease_7_29+sd), width=.2,
                position=position_dodge(.9))+theme(
                  axis.text.x = element_text(color = "black", size = 16),
                  axis.text.y = element_text(color = "black", size = 16),
                  axis.title.x = element_text(color = "black", size = 16),
                  axis.title.y = element_text(color = "black", size = 16)
                )

#respiration rates for inoculating cultures
# Microresp for each of the 5 weeks of inoculating (CommA, CommB, CommC)
#set as factor
orig_resp$week <- as.factor(orig_resp$week)
orig_resp$rep <- as.factor(orig_resp$rep)
orig_resp$treatment <- as.factor(orig_resp$treatment)

#calculate rate fo all the data
resp_ppm_rate <- orig_resp %>% mutate(rate24 = (orig_resp$T24 - orig_resp$T0) / 24)
ggplot(resp_ppm_rate, aes(x=treatment, y=rate24, color=week)) + geom_jitter()

#remove week 1 data, something not quite right with microresp plate
resp_ppm_rate_filt <- resp_ppm_rate %>% filter(!(week == 1))
resp_ppm_rate_filt <- resp_ppm_rate_filt[,-7]

#average reps for each week
resp_ppm_rate_filt_mean <- plyr::ddply(resp_ppm_rate_filt, c("treatment", "week"), summarise,
                                       N    = sum(!is.na(rate24)),
                                       ppm_mean = mean(rate24, na.rm = TRUE),
                                       sd   = sd(rate24, na.rm = TRUE),
                                       se   = sd / sqrt(N))


#average each week
resp_ave <- plyr::ddply(resp_ppm_rate_filt_mean, c("treatment"), summarise,
                        N    = sum(!is.na(ppm_mean)),
                        mean = mean(ppm_mean, na.rm = TRUE),
                        sd   = sd(ppm_mean, na.rm = TRUE),
                        se   = sd / sqrt(N))

glm.treatment_resp<- brm(rate24 ~ treatment, data= resp_ppm_rate_filt, family=Gamma(link="log"), iter=10000)
summary(glm.treatment_resp)
pp_check(glm.treatment_resp)

#posterior distributions from models for inoculating cultures
posterior_chit <- as.data.frame(glm.chit.org)
posterior_prot <- as.data.frame(glm.prot.org)
posterior_growth <- as.data.frame(glm.growth_cult)
posterior_resp <- as.data.frame(glm.treatment_resp)

pchit_pred <- ggpredict(glm.chit.org, "treatment")
pprot_pred <- ggpredict(glm.prot.org, "treatment")
pgrowth_pred <- ggpredict(glm.growth_cult, "treatment")
presp_pred <- ggpredict(glm.treatment_resp, "treatment")


pchit_melt <- posterior_chit[,1:3]
pprot_melt <- posterior_prot[,1:3]
pgrowth_melt <- posterior_growth[,1:3]
presp_melt <- posterior_resp[,1:3]

pchit_melt <- reshape2::melt(pchit_melt)
pprot_melt <- reshape2::melt(pprot_melt)
pgrowth_melt <- reshape2::melt(pgrowth_melt)
presp_melt <- reshape2::melt(presp_melt)

posterior1 <- mcmc_intervals_data(glm.chit.org, 
                                           prob_outer=0.95,
                                           prob=0.5)

posterior1$nonzero <- NA
posterior1$nonzero[posterior1$ll>0 & posterior1$hh>0] <- "nonzero"
posterior1$nonzero[posterior1$ll<0 & posterior1$hh<0] <- "nonzero"
posterior1$nonzero[is.na(posterior1$nonzero)] <- "zero"

level_order <- c("b_treatmentCommC", "b_treatmentCommB","b_Intercept") 

#FIGS2A
posterior1$parameter <- factor(posterior1$parameter, level = level_order)

posterior1<- posterior1[1:3,]
ggplot(posterior1, aes(x = parameter,
                       shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "black") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_color_manual(name="",
                     values = c("grey60", "#484c8d")) +
  scale_shape_manual(values=c(16, 17), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +theme_bw() + 
  theme(axis.text.y = element_text( size=7), 
        axis.text.x=element_text(size=7),
        axis.title = element_text(size=7), 
        legend.text = element_text(size=7)) +
  xlab(NULL) +
  ylab("Probability of effect of\nTreatment on Chitinase Activity")+
  guides(linetype=FALSE) +theme_classic()+theme(legend.position="none")

#prot orig cultures
posterior2 <- mcmc_intervals_data(glm.prot.org, 
                                  prob_outer=0.95,
                                  prob=0.5)

posterior2$nonzero <- NA
posterior2$nonzero[posterior2$ll>0 & posterior2$hh>0] <- "nonzero"
posterior2$nonzero[posterior2$ll<0 & posterior2$hh<0] <- "nonzero"
posterior2$nonzero[is.na(posterior2$nonzero)] <- "zero"

level_order <- c("b_treatmentCommC", "b_treatmentCommB","b_Intercept") 

posterior2$parameter <- factor(posterior2$parameter, level = level_order)

posterior2<- posterior2[1:3,]

#FIGS2B
ggplot(posterior2, aes(x = parameter,
                       shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "black") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_color_manual(name="",
                     values = c("grey60", "#484c8d")) +
  scale_shape_manual(values=c(16, 17), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +theme_bw() + 
  theme(axis.text.y = element_text( size=7), 
        axis.text.x=element_text(size=7),
        axis.title = element_text(size=7), 
        legend.text = element_text(size=7)) +
  xlab(NULL) +
  ylab("Probability of effect of\nTreatment on Protease Activity")+
  guides(linetype=FALSE)  +theme_classic()+theme(legend.position="none")

#resp orig cultures
posterior3 <- mcmc_intervals_data(glm.treatment_resp, 
                                  prob_outer=0.95,
                                  prob=0.5)

posterior3$nonzero <- NA
posterior3$nonzero[posterior3$ll>0 & posterior3$hh>0] <- "nonzero"
posterior3$nonzero[posterior3$ll<0 & posterior3$hh<0] <- "nonzero"
posterior3$nonzero[is.na(posterior3$nonzero)] <- "zero"

level_order <- c("b_treatmentCommC", "b_treatmentCommB","b_Intercept") 

posterior3$parameter <- factor(posterior3$parameter, level = level_order)

posterior3<- posterior3[1:3,]

#FIGS2D
ggplot(posterior3, aes(x = parameter,
                       shape=nonzero)) + 
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "black") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_color_manual(name="",
                     values = c("grey60", "#484c8d")) +
  scale_shape_manual(values=c(16, 17), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +theme_bw() + 
  theme(axis.text.y = element_text( size=7), 
        axis.text.x=element_text(size=7),
        axis.title = element_text(size=7), 
        legend.text = element_text(size=7)) +
  xlab(NULL) +
  ylab("Probability of effect of\nTreatment on Bacterial Respiration")+
  guides(linetype=FALSE) +theme_classic()+theme(legend.position="none")

#growth orig cultures
posterior4 <- mcmc_intervals_data(glm.growth_cult, 
                                  prob_outer=0.95,
                                  prob=0.5)

posterior4$nonzero <- NA
posterior4$nonzero[posterior4$ll>0 & posterior4$hh>0] <- "nonzero"
posterior4$nonzero[posterior4$ll<0 & posterior4$hh<0] <- "nonzero"
posterior4$nonzero[is.na(posterior4$nonzero)] <- "zero"

level_order <- c("b_treatmentCommC", "b_treatmentCommB","b_Intercept") 

posterior4$parameter <- factor(posterior4$parameter, level = level_order)

posterior4<- posterior4[1:3,]

#FIGS2C
ggplot(posterior4, aes(x = parameter,
                       shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "black") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_color_manual(name="",
                     values = c("grey60", "#484c8d")) +
  scale_shape_manual(values=c(16, 17), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +theme_bw() + 
  theme(axis.text.y = element_text( size=7), 
        axis.text.x=element_text(size=7),
        axis.title = element_text(size=7), 
        legend.text = element_text(size=7)) +
  xlab(NULL) +
  ylab("Probability of effect of\nTreatment on Bacterial Growth Rate")+
  guides(linetype=FALSE)+ theme_classic()+theme(legend.position="none")

############################ Physiological Data for all plants day 1 through 55 ############################
#### CHITINASE ACTIVITY ####
data_filt2 <- data_filt
data_filt2$chitinase_rate_uM_per_min <- ifelse(data_filt2$chitinase_rate_uM_per_min <= 0, 0.001, data_filt2$chitinase_rate_uM_per_min)
glm_treatment_chitall<- brm(chitinase_rate_uM_per_min ~ treatment+ (1|day), data= data_filt2,
family = Gamma(link = "log"), iter = 10000, chains = 4, cores = 4,control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth=20))
saveRDS(glm_treatment_chitall, "glm_treatment_chitall.RDS")
glm_treatment_chitall<-readRDS("glm_treatment_chitall.RDS")
pp_check(glm_treatment_chitall)

#Subset to only those samples in metatranscriptomic analysis
col_treat <- c("#014b7a", "#ffaf49", "#44b7c2")

## ALL TIME POINTS, ALL METATRANSCRIPTOMIC PLANTS ONLY##
tran_chit <- data_filt %>%  subset(., plant_number %in% c("5", "18", "45", "19", "43", "1", "13", "3", "33"))
#glm_tran_chit<- brm(chitinase_rate_uM_per_min ~ treatment + (1|day), data= tran_chit, family=Gamma(link="log"), iter = 10000, chains = 4, cores = 4)
#saveRDS(glm_tran_chit, "glm_tran_chit.RDS")
glm_tran_chit<-readRDS("glm_tran_chit.RDS")

summary(glm_tran_chit)
posterior_interval(glm_tran_chit)
pmchit <- ggeffects::ggpredict(glm_tran_chit, "treatment")
colnames(pmchit)[1] <- "treatment"
colnames(pmchit)[2] <- "chitinase_rate_uM_per_min"
level_order2 <- c("CommA", "CommB", "CommC") 

##FIGURE 6A
ggplot(tran_chit, 
       mapping = aes(x = factor(treatment, level=level_order2), 
                     y = (chitinase_rate_uM_per_min*1000), 
                     fill = treatment, color=treatment, group=treatment))  +
  geom_boxplot(outlier.shape=NA, alpha=.5)+ geom_jitter(width=.4, size=3)+
  scale_color_manual(values=c("#004488", "#ffaf49", "#44b7c2")) +
  scale_fill_manual(values=c("#004488", "#ffaf49", "#44b7c2")) +
  geom_point(data=pmchit, aes(x=treatment, y=(chitinase_rate_uM_per_min*1000)), color="black", size=5) +
  geom_linerange(data=pmchit,aes(ymin=(conf.low*1000), ymax=(conf.high*1000)), size=2, color="black",
                 position=position_dodge(width = 0.5)) + labs(
                   x = "Treatment", 
                   y = "Chitinase Activity") + theme_classic()+
  theme(legend.position="none")+theme(
    axis.text.x = element_text(color = "black", size = 24),
    axis.text.y = element_text(color = "black", size = 24),
    axis.title.x = element_text(color = "black", size = 24),
    axis.title.y = element_text(color = "black", size = 24)
  )
posterior1 <- mcmc_intervals_data(glm_tran_chit, 
                                  prob_outer=0.95,
                                  prob=0.5)

posterior1$nonzero <- NA
posterior1$nonzero[posterior1$ll>0 & posterior1$hh>0] <- "nonzero"
posterior1$nonzero[posterior1$ll<0 & posterior1$hh<0] <- "nonzero"
posterior1$nonzero[is.na(posterior1$nonzero)] <- "zero"

level_order <- c("b_treatmentCommC", "b_treatmentCommB","b_Intercept") 

#FIGS6A
posterior1$parameter <- factor(posterior1$parameter, level = level_order)

posterior1<- posterior1[1:3,]
ggplot(posterior1, aes(x = parameter,
                       shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "black") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_color_manual(name="",
                     values = c("grey60", "#484c8d")) +
  scale_shape_manual(values=c(16, 17), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +theme_bw() + 
  theme(axis.text.y = element_text( size=7), 
        axis.text.x=element_text(size=7),
        axis.title = element_text(size=7), 
        legend.text = element_text(size=7)) +
  xlab(NULL) +
  ylab("Probability of effect of\nTreatment on Chitinase Activity")+
  guides(linetype=FALSE) +theme_classic()+theme(legend.position="none")

#### PROTEASE ACTIVITY ####
## ALL TIME POINTS, ALL METATRANSCRIPTOMIC PLANTS ONLY##
tran_prot <- data_filt %>%  subset(., plant_number %in% c("5", "18", "45", "19", "43", "1", "13", "3", "33"))
#glm_tran_prot<- brm(protease_rate_nM_per_min ~treatment +(1|day), data= tran_prot, family=Gamma(link="log"), iter = 10000, chains = 4, cores = 4)
#saveRDS(glm_tran_prot, "glm_tran_prot.RDS")
glm_tran_prot<-readRDS("glm_tran_prot.RDS")

summary(glm_tran_prot)

pmprot <- ggeffects::ggpredict(glm_tran_prot, "treatment")
colnames(pmprot)[1] <- "treatment"
colnames(pmprot)[2] <- "protease_rate_nM_per_min"

#Fig6B
ggplot(tran_prot, 
       mapping = aes(x = factor(treatment, level=level_order2), 
                     y = protease_rate_nM_per_min, 
                     fill = treatment, color=treatment))  +
  geom_boxplot(outlier.shape=NA, alpha=.5) + geom_jitter(width=.4, size=3)+
  scale_color_manual(values=c("#004488", "#ffaf49", "#44b7c2")) +
  scale_fill_manual(values=c("#004488", "#ffaf49", "#44b7c2")) +
  geom_point(data=pmprot, aes(x=treatment, y=protease_rate_nM_per_min), color="black", size=5) +
  geom_linerange(data=pmprot,aes(ymin=conf.low, ymax=conf.high), size=2, color="black",
                 position=position_dodge(width = 0.5)) + labs(
                   x = "Treatment", 
                   y = "Protease Activity") + theme_classic()+
  theme(legend.position="none")+theme(
    axis.text.x = element_text(color = "black", size = 24),
    axis.text.y = element_text(color = "black", size = 24),
    axis.title.x = element_text(color = "black", size = 24),
    axis.title.y = element_text(color = "black", size = 24)
  )

posterior1 <- mcmc_intervals_data(glm_tran_prot, 
                                  prob_outer=0.95,
                                  prob=0.5)

posterior1$nonzero <- NA
posterior1$nonzero[posterior1$ll>0 & posterior1$hh>0] <- "nonzero"
posterior1$nonzero[posterior1$ll<0 & posterior1$hh<0] <- "nonzero"
posterior1$nonzero[is.na(posterior1$nonzero)] <- "zero"

level_order <- c("b_treatmentCommC", "b_treatmentCommB","b_Intercept") 

#FIGS6B
posterior1$parameter <- factor(posterior1$parameter, level = level_order)

posterior1<- posterior1[1:3,]
ggplot(posterior1, aes(x = parameter,
                       shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "black") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_color_manual(name="",
                     values = c("grey60", "#484c8d")) +
  scale_shape_manual(values=c(16, 17), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +theme_bw() + 
  theme(axis.text.y = element_text( size=7), 
        axis.text.x=element_text(size=7),
        axis.title = element_text(size=7), 
        legend.text = element_text(size=7)) +
  xlab(NULL) +
  ylab("Probability of effect of\nTreatment on Protease Activity")+
  guides(linetype=FALSE) +theme_classic()+theme(legend.position="none")


#### RESPIRATION ####
#glm_treatment_respall<- brm(bacterial_respiration_rate_at_24hr ~ treatment+ (1|day), data= data_filt2,
#family = Gamma(link = "log"), iter = 10000, chains = 4, cores = 4,control = list(adapt_delta = 0.999, stepsize = 0.001, max_treedepth=20))
#saveRDS(glm_treatment_respall, "glm_treatment_respall.RDS")
glm_treatment_respall<-readRDS("glm_treatment_respall.RDS")
mcmc_plot(glm_treatment_respall)

#### GROWTH RATE ####
glm_treatment_growth<- brm(wk8_r ~ treatment, data= data_filt, family=Gamma(link="log"), iter = 20000, chains = 4, cores = 4, control = list(adapt_delta = 0.999, max_treedepth = 20))
summary(glm_treatment_growth)
mcmc_plot(glm_treatment_growth)

#### PLANT TRAITS ####
data_filt$nitrogen1_content_grams <- (data_filt$leaf1_percent_N / 100) * data_filt$drymass_leaf1_grams
data_filt$nitrogen2_content_grams <- (data_filt$leaf2_percent_N / 100) * data_filt$drymass_leaf2_grams
data_filt$carbon1_content_grams <- (data_filt$leaf1_percent_C / 100) * data_filt$drymass_leaf1_grams
data_filt$carbon2_content_grams <- (data_filt$leaf2_percent_C / 100) * data_filt$drymass_leaf2_grams

detach("package:plyr", unload=TRUE)
filtered_data <- data_filt %>% 
  filter(!is.na(nitrogen2_content_grams))
cor_by_treatment <- filtered_data %>%
  dplyr::group_by(treatment) %>%
  summarize(correlation_squared = cor(nitrogen1_content_grams, nitrogen2_content_grams)^2)

# pitcher carbon content
#FigS4C
ggplot(data=data_filt,aes(x = treatment, y = carbon1_content_grams, group_by=treatment,fill = treatment)) +geom_boxplot(outlier.shape=NA, size=1, alpha=.7)+
  geom_jitter(aes(color=treatment),size=3) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_color_manual(values=c("#5C4033", "#014b7a", "#ffaf49", "#44b7c2", "gray"))+
  scale_fill_manual(values=c("#5C4033", "#014b7a", "#ffaf49", "#44b7c2", "gray"))
level_order2 <- c("CommA", "CommB", "CommC", "ACM", "WATER") 

# pitcher length
day1_len <- filter(data_filt, day == "1")
day8_len <- filter(data_filt, day == "8")
day15_len <- filter(data_filt, day == "15")

glm_treatment_lengthday1<- brm(pitcher_length_cm ~ treatment, data= day1_len, family=Gamma(link="log"), iter = 20000, chains = 4, cores = 4, control = list(adapt_delta = 0.999, max_treedepth = 20))
glm_treatment_lengthday8<- brm(pitcher_length_cm ~ treatment, data= day8_len, family=Gamma(link="log"), iter = 20000, chains = 4, cores = 4, control = list(adapt_delta = 0.999, max_treedepth = 20))
glm_treatment_lengthday15<- brm(pitcher_length_cm ~ treatment, data= day15_len, family=Gamma(link="log"), iter = 20000, chains = 4, cores = 4, control = list(adapt_delta = 0.999, max_treedepth = 20))

summary(glm_treatment_lengthday1)
mcmc_plot(glm_treatment_lengthday1, regex_pars="b_", prob=.95, prob_outer=.95)+
  theme_classic() 
mcmc_plot(glm_treatment_lengthday8, regex_pars="b_", prob=.95, prob_outer=.95)+
  theme_classic() 
mcmc_plot(glm_treatment_lengthday15, regex_pars="b_", prob=.95, prob_outer=.95)+
  theme_classic() 



glm_treatment_length2<- brm(pitcher_length_cm ~ treatment+week+(1|plant_number), data= data_filt, family=Gamma(link="log"), iter = 20000, chains = 4, cores = 4, control = list(adapt_delta = 0.999, max_treedepth = 20))
mcmc_plot(glm_treatment_length2)

pm1a <- ggeffects::ggpredict(glm_treatment_length2, c("week", "treatment"))
colnames(pm1a)[1] <- "week"
colnames(pm1a)[2] <- "pitcher_length_cm"
colnames(pm1a)[5] <- "treatment"

ggplot() +
  geom_jitter(data = data_filt, aes(x = week, y = pitcher_length_cm, color = factor(treatment, level=level_order2), fill = factor(treatment, level=level_order), width=.1)) +
  geom_line(data = pm1a, aes(x = week, y = pitcher_length_cm, color = treatment, group = treatment), lwd=1.5) +
  geom_ribbon(data = pm1a, aes(x = week, ymin = conf.low, ymax = conf.high, fill = treatment), alpha = 0.1) + # Shaded area for confidence interval
  scale_color_manual(values = c("#014b7a", "#ffaf49", "#44b7c2", "#5C4033", "gray")) +
  scale_fill_manual(values = c( "#014b7a", "#ffaf49", "#44b7c2", "#5C4033","gray"))+theme_classic()+
  theme(legend.position ="none")+ theme(text = element_text(color = "black", size = 16 ))+
  xlab("Week")+ ylab("Pitcher Length (cm)")

posterior_len <- as.data.frame(glm_treatment_length2)
posterior_len_melt <- posterior_len[,1:5]
posterior_len_melt <- reshape2::melt(posterior_len_melt)
level_order <- c("b_treatmentWATER", "b_treatmentCommC", "b_treatmentCommB", "b_treatmentCommA", "b_Intercept") 

ggplot(posterior_len_melt, aes(x = value, y =factor(variable, level = level_order),
                                    fill = variable)) +
  stat_halfeye(slab_alpha=0.50, .width=c(.95, 0.95),alpha=.6) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + xlab("")+
  scale_fill_manual(values=c("#5C4033","#004488", "#ffaf49", "#44b7c2", "gray"))

# pitcher photosynthetic rate
#FigS7A
ggplot(data=data_filt,aes(x = treatment, y = photosynthetic_rate_leaf1_μmol_per_m2sec, group_by=treatment,fill = treatment)) +geom_boxplot(outlier.shape=NA, size=1, alpha=.7)+
  geom_jitter(aes(color=treatment),size=3) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_color_manual(values=c("gray","#5C4033", "#014b7a", "#ffaf49", "#44b7c2"))+
  scale_fill_manual(values=c("gray","#5C4033", "#014b7a", "#ffaf49", "#44b7c2"))

data_filt$photosynthetic_rate_leaf1_μmol_per_m2sec <- ifelse(data_filt$photosynthetic_rate_leaf1_μmol_per_m2sec <= 0, 0.0001, data_filt$photosynthetic_rate_leaf1_μmol_per_m2sec)

glm_treatment_leaf1_photo<- brm(photosynthetic_rate_leaf1_μmol_per_m2sec ~ treatment, data= data_filt, family=Gamma(link="log"), iter = 20000, chains = 4, cores = 4, control = list(adapt_delta = 0.999, max_treedepth = 20))
summary(glm_treatment_leaf1_photo)
mcmc_plot(glm_treatment_leaf1_photo)
posterior6 <- mcmc_intervals_data(glm_treatment_leaf1_photo, 
                                  prob_outer=0.95,
                                  prob=0.5)

posterior6$nonzero <- NA
posterior6$nonzero[posterior6$ll>0 & posterior6$hh>0] <- "nonzero"
posterior6$nonzero[posterior6$ll<0 & posterior6$hh<0] <- "nonzero"
posterior6$nonzero[is.na(posterior6$nonzero)] <- "zero"
posterior6<- posterior6[1:5,]
level_order <- c("b_treatmentWATER", "b_treatmentCommC", "b_treatmentCommB","b_treatmentCommA", "b_Intercept") 
posterior6$parameter <- factor(posterior6$parameter, level = level_order)

#FigS7B
ggplot(posterior6, aes(x = parameter,
                       shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_color_manual(name="",
                     values = c("grey60", "#484c8d")) +
  scale_shape_manual(values=c(16, 17), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +theme_bw() + 
  theme(axis.text.y = element_text( size=7), 
        axis.text.x=element_text(size=7),
        axis.title = element_text(size=7), 
        legend.text = element_text(size=7)) +
  xlab(NULL) +
  ylab("Probability of effect of\non Pitcher Photosynthesis ")+
  guides(linetype=FALSE) + theme(legend.position = "none")

# pitcher net leaf production
#FigS7G
ggplot(data=data_filt,aes(x = treatment, y = new_leaves, group_by=treatment,fill = treatment)) +geom_boxplot(outlier.shape=NA, size=1, alpha=.7)+
  geom_jitter(aes(color=treatment),size=3) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_color_manual(values=c("#5C4033", "#014b7a", "#ffaf49", "#44b7c2", "gray"))+
  scale_fill_manual(values=c("#5C4033", "#014b7a", "#ffaf49", "#44b7c2", "gray"))


glm_treatment_newleaves<- brm(new_leaves ~ treatment, data= data_filt)
pleaves <- ggeffects::ggpredict(glm_treatment_newleaves, terms = c("treatment"))

summary(glm_treatment_newleaves)
posterior7 <- mcmc_intervals_data(glm_treatment_newleaves, 
                                  prob_outer=0.95,
                                  prob=0.5)

posterior7$nonzero <- NA
posterior7$nonzero[posterior7$ll>0 & posterior7$hh>0] <- "nonzero"
posterior7$nonzero[posterior7$ll<0 & posterior7$hh<0] <- "nonzero"
posterior7$nonzero[is.na(posterior7$nonzero)] <- "zero"
posterior7<- posterior7[1:5,]
level_order <- c("b_treatmentWATER", "b_treatmentCommC", "b_treatmentCommB","b_treatmentCommA", "b_Intercept") 
posterior7$parameter <- factor(posterior7$parameter, level = level_order)

#percent positive for CommC
posterior_glm_treatment_newleaves <- as.data.frame(glm_treatment_newleaves)
CommC_leaves <- posterior_glm_treatment_newleaves %>% filter(b_treatmentCommC <0)
nrow(CommC_leaves)/nrow(posterior_glm_treatment_newleaves) #0.9655

#estimate of effect
glm_treatment_newleaves
exp(1.40)

#FigS7H
ggplot(posterior7, aes(x = parameter,
                       shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_color_manual(name="",
                     values = c("grey60", "#484c8d")) +
  scale_shape_manual(values=c(17), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +theme_bw() + 
  theme(axis.text.y = element_text( size=7), 
        axis.text.x=element_text(size=7),
        axis.title = element_text(size=7), 
        legend.text = element_text(size=7)) +
  xlab(NULL) +
  ylab("Estimated effect on pitchers gained/lost")+
  guides(linetype=FALSE) + theme(legend.position = "none")


#### Correlation between plant traits ####
#subset data to just plant traits
plant_trait <- data.frame(data_filt)
plant_trait <- plant_trait[,c(13,14,17,18,19,23,24,25,26,29,30,31,32,39,40,41,53,54,55,56)]

plant_trait_naomit <- na.omit(plant_trait)
data_normalized <- scale(plant_trait_naomit)
corr_matrix <- cor(data_normalized)
ggcorrplot(corr_matrix)
data.pca <- princomp(corr_matrix)
summary(data.pca)
data.pca$loadings[, 1:2]
fviz_eig(data.pca, addlabels = TRUE)
fviz_pca_var(data.pca, col.var = "black")
fviz_cos2(data.pca, choice = "var", axes = 1:2)

#FigS9
fviz_pca_var(data.pca, col.var = "x",
             gradient.cols = c("navy","black", "maroon"),
             repel = TRUE)


#### ECOPLATES ####
# no blanks, experimental controls, or damaged plants
#removed plants because after initial analysis these were extreme outliers
#"21", "28", "41" CommC
#"15", "34" CommA
#"27" CommB
eco_filt <- subset(ecodat, treatment != "WATER" & treatment != "ACM" & 
                 !((day == "55") & (plant %in% c("25", "44", "51", "9", "29", "21", "41", "27", "34", "28","15"))))


com5 <- eco_filt

com5_meta <- subset(com5, select= c(plant, treatment, day))
com5 <- subset(com5, select= -c(plant, treatment, day))
com5[com5 < 0] <- 0

#run NMDS on ecoplate data, day 1 and day 55
set.seed(123)
nmds5 <- metaMDS(com5, k=2, trymax=500)

#get the data scores from the nmds5
all.scores <- as.data.frame(scores(nmds5, "sites"))
#add metadata to the scores info
all.scores <- cbind(all.scores, com5_meta)
all.scores$day <- as.factor(all.scores$day)

# make a new column that has unique names for week and treatment
all.scores <- all.scores %>% 
  unite(hull.id, treatment, day, remove = FALSE)

hullsall <- all.scores %>%
  group_by(hull.id) %>%
  slice(chull(NMDS1, NMDS2))

#FigS3
ggplot(all.scores, aes(x = NMDS1, y = NMDS2, label=plant)) + 
  geom_point(size = 4, aes(colour = treatment, shape = day)) +
  theme(axis.text.y = element_text(colour = "black", size = 16, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 16), 
        legend.text = element_text(size = 16, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 16), 
        axis.title.x = element_text(face = "bold", size = 16, colour = "black"), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "treatment", y = "NMDS2", shape = "day")  + 
  scale_colour_manual(values = c("#014b7a", "#ffaf49", "#44b7c2"))  + aes(fill = factor(hull.id)) + geom_polygon(data = hullsall, alpha = 0.2) + guides(fill=FALSE)+ 
  scale_fill_manual(values = c("#014b7a","#014b7a", "#ffaf49", "#ffaf49", "#44b7c2","#44b7c2"))+geom_text()

#day 1 only
com1 <- subset(ecodat, treatment != "WATER" & treatment != "ACM" & day != "1")
com1_meta <- subset(com1, select= c(plant, treatment, day))
com1 <- subset(com1, select= -c(plant, treatment, day))
com1[com1 < 0] <- 0

#run NMDS on day 1 data only
nmds1 <- metaMDS(com1, k=2, trymax=500)
ordiplot(nmds1, display = "sites", type = "t")
ordiplot(nmds1, display = "species", type = "t")

wk1.scores <- as.data.frame(scores(nmds1, "sites"))
wk1.scores <- cbind(wk1.scores, com1_meta)
wk1.scores$day <- as.factor(wk1.scores$day)

# make a new column that has unique names for week and treatment
wk1.scores <- wk1.scores %>% 
  unite(hull.id, treatment, day, remove = FALSE)

hullswk1 <- wk1.scores %>%
  group_by(hull.id) %>%
  slice(chull(NMDS1, NMDS2))



###add it prior inoculant samples one each
#day 1 only
com1 <- subset(ecodat, treatment != "WATER" & treatment != "ACM" & day != "55")
com2 <- read.csv("EcoPlates_06_17_2021.csv")
com_all <- rbind(com1, com2)
com_all_meta <- subset(com_all, select= c(plant, treatment, day))
com_all <- subset(com_all, select= -c(plant, treatment, day))
com_all[com_all < 0] <- 0

#run NMDS on day 1 data only
nmds1 <- metaMDS(com_all, k=2, trymax=500)#0.1319898

wk1.scores <- as.data.frame(scores(nmds1, "sites"))
wk1.scores <- cbind(wk1.scores, com_all_meta)
wk1.scores$day <- as.factor(wk1.scores$day)

# make a new column that has unique names for week and treatment
wk1.scores <- wk1.scores %>% 
  unite(hull.id, treatment, day, remove = FALSE)

hullswk1 <- wk1.scores %>%
  group_by(hull.id) %>%
  slice(chull(NMDS1, NMDS2))

#Fig1C
ggplot(wk1.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = treatment, shape = day)) +
  theme(axis.text.y = element_text(colour = "black", size = 16), 
        axis.text.x = element_text(colour = "black",  size = 16), 
        legend.text = element_text(size = 16,  colour ="black"), 
        legend.position = "right", axis.title.y = element_text( size = 16), 
        axis.title.x = element_text(size = 16, colour = "black"), 
        legend.title = element_text(size = 16, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "treatment", y = "NMDS2", shape = "day")  + 
  scale_colour_manual(values = c("#014b7a", "#ffaf49", "#44b7c2", "#014b7a", "#ffaf49", "#44b7c2"))  + aes(fill = factor(hull.id)) + geom_polygon(data = hullswk1, alpha = 0.2) + guides(fill=FALSE)+ 
  scale_fill_manual(values = c("#014b7a", "#ffaf49", "#44b7c2", "#014b7a", "#ffaf49", "#44b7c2"))

com1_meta <- subset(com1, select= c(plant, treatment, day))
com1 <- subset(com1, select= -c(plant, treatment, day))
com1[com1 < 0] <- 0

## Ecoplate Statistical Analysis
# beta disper to check dispersion day 1 and day 55
ecoall.bd <- betadisper(vegdist(com5), com5_meta$treatment)
ecoall.bd
boxplot(ecoall.bd)
anova(ecoall.bd)# p=0.5602

ecoall.bdday <- betadisper(vegdist(com5), com5_meta$day)
ecoall.bdday
boxplot(ecoall.bdday)
anova(ecoall.bdday)# p=0.06473

# beta disper to check dispersion day 8 and day 55
ecoall.day1 <- betadisper(vegdist(com1), com1_meta$treatment)
ecoall.day1
boxplot(ecoall.day1)
anova(ecoall.day1)# p=0.02602

# PERMANOVA for categorical variables (factors)
#week1
adonis.eco.w1 <- adonis2(com1 ~ treatment, data=wk1.scores, method = "bray")
adonis.eco.w1

#day 8 and day 55
# set number of permutations
perm <- how(nperm = 999)
#specify a random variable (plant).
setBlocks(perm) <- with(all.scores, plant)
adonis <- adonis2(com5 ~ treatment * day, data=all.scores, permutations = perm, method = "bray")
adonis

# pairwise adonis
ecoall.ad.pw <- pairwise.adonis2(com5 ~ treatment * day, data = com5_meta, strata = 'plant')
ecoall.ad.pw

#### ECOPLATES THROUGH TIME ####
#week one doesn't have plant 35 because that pitcher had low volume
dat.1 <- eco_filt[!(eco_filt$plant=="35"),]
dat.1 <- replace(dat.1, dat.1 < 0, 0.001)

dat.carbons <- melt(dat.1[, c("plant","treatment", "day", 
                              "Pyruvic.acid.methyl.ester","Tween.40","Tween.80","Alpha.cyclodextrin",
                              "Glycogen","D.Cellobiose","alpha.D.Lactose","Beta.Methyl.D.Glucoside",
                              "D.Xylose","I.Erythritol","D.Mannitol","N.Acetyl.D.Glucosamine",
                              "D.Glucosaminic.acid","Glucose.1.phosphate","D.L.alpha.glycerol.phosphate",
                              "D.Galactonic.acid.gamma.Lactone","D.Galacturonic.acid","X2.Hydroxy.benzoic.acid",
                              "X4.Hydroxy.benzoic.acid","gamma.Hydroxybutyric.acid","Itaconic.acid",
                              "alpha.Ketobutyric.acid","D.Malic.acid","L.Arginine","L.Asparagine",
                              "L.Phenylalanine","L.Serine","L.Threonine","Glycyl.L.glutamic.acid",
                              "Phenylethylamine","Putrescine")], 
                    measure.vars = c("Pyruvic.acid.methyl.ester","Tween.40","Tween.80","Alpha.cyclodextrin",
                                     "Glycogen","D.Cellobiose","alpha.D.Lactose","Beta.Methyl.D.Glucoside",
                                     "D.Xylose","I.Erythritol","D.Mannitol","N.Acetyl.D.Glucosamine",
                                     "D.Glucosaminic.acid","Glucose.1.phosphate","D.L.alpha.glycerol.phosphate",
                                     "D.Galactonic.acid.gamma.Lactone","D.Galacturonic.acid","X2.Hydroxy.benzoic.acid",
                                     "X4.Hydroxy.benzoic.acid","gamma.Hydroxybutyric.acid","Itaconic.acid",
                                     "alpha.Ketobutyric.acid","D.Malic.acid","L.Arginine","L.Asparagine",
                                     "L.Phenylalanine","L.Serine","L.Threonine","Glycyl.L.glutamic.acid",
                                     "Phenylethylamine","Putrescine"))
dat.carbons <- within(dat.carbons, {
  c2 <- as.integer(variable == "Pyruvic.acid.methyl.ester")
  c3<- as.integer(variable == "Tween.40")
  c4<- as.integer(variable == "Tween.80")
  c5<- as.integer(variable == "Alpha.cyclodextrin")
  c6<- as.integer(variable == "Glycogen")
  c7<- as.integer(variable == "D.Cellobiose")
  c8<- as.integer(variable == "alpha.D.Lactose")
  c9<- as.integer(variable == "Beta.Methyl.D.Glucoside")
  c10<- as.integer(variable == "D.Xylose")
  c11<- as.integer(variable == "I.Erythritol")
  c12<- as.integer(variable == "D.Mannitol")
  c13<- as.integer(variable == "N.Acetyl.D.Glucosamine")
  c14<- as.integer(variable == "D.Glucosaminic.acid")
  c15<- as.integer(variable == "Glucose.1.Phosphate")
  c16<- as.integer(variable == "D.L.alpha.glycerol.phosphate")
  c17<- as.integer(variable == "D.Galactonic.acid.gamma.Lactone")
  c18<- as.integer(variable == "D.Galacturonic.acid")
  c19<- as.integer(variable == "X2.Hydroxy.benzoic.acid")
  c20<- as.integer(variable == "X4.Hydroxy.benzoic.acid")
  c21<- as.integer(variable == "gamma.Hydroxybutyric.acid")
  c22<- as.integer(variable == "Itaconic.acid")
  c23<- as.integer(variable == "alpha.Ketobutyric.acid")
  c24<- as.integer(variable == "D.Malic.acid")
  c25<- as.integer(variable == "L.Arginine")
  c26<- as.integer(variable == "L.Asparagine")
  c27<- as.integer(variable == "L.Phenylalanine")
  c28<- as.integer(variable == "L.Serine")
  c29<- as.integer(variable == "L.Threonine")
  c30<- as.integer(variable == "Glycyl.L.glutamic.acid")
  c31<- as.integer(variable == "Phenylethylamine")
  c32<- as.integer(variable == "Putrescine")
})

#model of treatment*substrate on nitrogen
dat.carbons2 <- dat.carbons %>%
  dplyr::rename("plant_number" = "plant")
data_filteco <- data_filt[,c(4,6, 8)]
data_filteco2 <- subset(data_filteco, day %in% c("1", "55"))
#data_filteco <- subset(data_filteco, day == "55")
dat.carbons2$plant_number <- as.character(dat.carbons2$plant_number)
dat.carbons2$day <- as.factor(dat.carbons2$day)

eco_nit <- left_join(dat.carbons2, data_filteco2, by = c("plant_number", "day"))


eco_nit <- eco_nit[,c(1:5)]

unique_sub <- unique(eco_nit$variable)
commB_summary <- aggregate(value ~ variable, data = eco_nit[eco_nit$treatment == "CommB", ], FUN = mean)  
sorted_sub <- commB_summary[order(commB_summary$value, decreasing = TRUE), "variable"]

eco_nit$variable <- factor(eco_nit$variable, levels = sorted_sub)
eco_nit$variable <- factor(eco_nit$variable, levels = sorted_sub)
detach("package:plyr", unload=TRUE)

means <- eco_nit %>%
  group_by(variable, treatment.x) %>%
  summarise(mean_value = mean(value))

means$treatment.x <- factor(means$treatment.x, levels = c("CommC", "CommB", "CommA"))

eco_nit$variable <- relevel(eco_nit$variable, ref = "X2.Hydroxy.benzoic.acid")


#ecom3 <- brm(value~treatment.x*variable, data=eco_nit , family=Gamma(link="log"), iter = 10000, chains = 4, cores = 4)
#saveRDS(ecom3, "Ecoplate_treatment.RDS")

ecom3 <- readRDS("Ecoplate_treatment.RDS")
mcmc_areas(ecom3, 
           prob = 0.90,
           prob_outer = 0.90, 
           point_est = "mean")+ theme_classic() + xaxis_text(color="black")+ yaxis_text(color="black") +
  geom_vline(xintercept=0, linetype="dotted", color="darkgray")

posterioreco <- mcmc_intervals_data(ecom3, 
                                  prob_outer=0.90,
                                  prob=0.5)

posterioreco$nonzero <- NA
posterioreco$nonzero[posterioreco$ll>0 & posterioreco$hh>0] <- "nonzero"
posterioreco$nonzero[posterioreco$ll<0 & posterioreco$hh<0] <- "nonzero"
posterioreco$nonzero[is.na(posterioreco$nonzero)] <- "zero"
posterioreco<- posterioreco[34:93,]
ggplot(posterioreco, aes(x = parameter,shape=nonzero)) + geom_hline(yintercept = 0, linetype = 3, 
size=1, color = "#b0b5b3") +  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
position= position_dodge(width=0.75),size = 3/4) +scale_color_manual(name="",values = c("grey60", "#484c8d")) +
scale_shape_manual(values=c(16, 17), labels=c("95% CI does\nnot contain zero", "95% CI\ncontains zero"))+
coord_flip() +theme_bw() + theme(axis.text.y = element_text( size=7), axis.text.x=element_text(size=7),
 axis.title = element_text(size=7), legend.text = element_text(size=7)) +xlab(NULL) +
  guides(linetype=FALSE)

treat_col <- c("#004488", "#ffaf49", "#44b7c2")

pecom3 <- ggpredict(ecom3, terms=c("treatment.x", "variable"))
colnames(pecom3)[1] <- "treatment.x"
colnames(pecom3)[2] <- "value"
colnames(pecom3)[6] <- "variable"

#Filter down to just those that look different from eachother
eco_nit2 <- filter(eco_nit, variable %in% c("N.Acetyl.D.Glucosamine", "Beta.Methyl.D.Glucoside", "D.Cellobiose",
                                  "Pyruvic.acid.methyl.ester", "Tween.80", "D.Galactonic.acid.gamma.Lactone",
                                  "L.Serine", "L.Threonine", "D.Glucosaminic.acid"))
pecom32 <- filter(pecom3, variable %in% c("N.Acetyl.D.Glucosamine", "Beta.Methyl.D.Glucoside", "D.Cellobiose",
                                            "Pyruvic.acid.methyl.ester", "Tween.80", "D.Galactonic.acid.gamma.Lactone",
                                            "L.Serine", "L.Threonine", "D.Glucosaminic.acid"))
eco_nit2$variable <- factor(eco_nit2$variable, levels = c("N.Acetyl.D.Glucosamine", "Beta.Methyl.D.Glucoside", "D.Cellobiose",
                                                          "Pyruvic.acid.methyl.ester", "Tween.80", "L.Serine", "L.Threonine","D.Galactonic.acid.gamma.Lactone",
                                                          "D.Glucosaminic.acid"))
ggplot(eco_nit2, 
       mapping = aes(x = variable, 
                     y = value, 
                     fill = treatment.x, color = treatment.x)) +
  geom_boxplot(outlier.shape = NA, alpha = .5, position = position_dodge(width = .9)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = .9), size = 1) +
  scale_color_manual(values = c("#004488", "#ffaf49", "#44b7c2")) +
  scale_fill_manual(values = c("#004488", "#ffaf49", "#44b7c2")) +
  geom_point(data = pecom32, aes(x = variable, y = value), color = "black", size = 1, 
             position = position_dodge(width = .9)) +
  geom_linerange(data = pecom32, aes(ymin = conf.low, ymax = conf.high), size = 1, color = "black",
                 position = position_dodge(width = .9)) + 
  theme_classic() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## VISUALIZATIONS
# Plot average value and se for each substrate within each treatment
# x-axis is the day1 timepoint and y-axis is the day55 timepoint.

#day 1 data
dat.carbons1<-dat.carbons[!(dat.carbons$day=="55"),]
avg.subst1 <- with(dat.carbons1 , aggregate(value, list("Treatment"=as.factor(treatment),"Substrate"=as.factor(variable)), 
                                            mean))
avg.subst1$se <- with(dat.carbons1, aggregate(value, list("Treatment"=as.factor(treatment),"Substrate"=as.factor(variable)),
                                              function(x) std.error(x)))[,3]

#day 55 data
dat.carbons2<-dat.carbons[!(dat.carbons$day=="1"),]
avg.subst2 <- with(dat.carbons2 , aggregate(value, list("Treatment"=as.factor(treatment),"Substrate"=as.factor(variable)), 
                                            mean))
avg.subst2$se <- with(dat.carbons2, aggregate(value, list("Treatment"=as.factor(treatment),"Substrate"=as.factor(variable)),
                                              function(x) std.error(x)))[,3]
avg.subst1$day8 <-avg.subst2$x
avg.subst1$day8.se <- avg.subst2$se

agg <- plyr::ddply(dat.carbons, c("treatment", "day", "variable"), summarise,
             N    = sum(!is.na(value)),
             mean = mean(value, na.rm = TRUE),
             sd   = sd(value, na.rm = TRUE),
             se   = sd / sqrt(N))

pd <- position_dodge(0.4)

pairs <- c(0:30, 0:30, 31:61, 31:61, 62:92, 62:92)

agg <- cbind(agg, pairs)

#lets look at the treatments individually
treatCommA <- subset(agg, treatment == "CommA")
treatCommB <- subset(agg, treatment == "CommB")
treatCommC <- subset(agg, treatment == "CommC")

### correlations
wide <- cbind(treatCommA, treatCommB)
wide <- cbind(wide, treatCommC)

# Define new column names
new_colnames <- c(
  "CommA", "CommA_day", "CommA_variable", "CommA_N", "CommA_mean", "CommA_sd", "CommA_se", "CommA_pairs",
  "CommB", "CommB_day", "CommB_variable", "CommB_N", "CommB_mean", "CommB_sd", "CommB_se", "CommB_pairs",
  "CommC", "CommC_day", "CommC_variable", "CommC_N", "CommC_mean", "CommC_sd", "CommC_se", "CommC_pairs"
)

# Assign the new column names
colnames(wide) <- new_colnames

wide$CommA_day <- as.factor(wide$CommA_day)
wide$CommB_day <- as.factor(wide$CommB_day)
wide$CommC_day <- as.factor(wide$CommC_day)

##CommBvsCommA
ggplot(data=wide, aes(x=CommA_mean, y=CommB_mean, color = CommA_day, label=CommA_variable))+geom_point()+
  geom_smooth(method="lm",aes(x=CommA_mean,y=CommB_mean))+geom_errorbar(aes(ymax=CommB_mean+CommB_se, ymin=CommB_mean-CommB_se,x=CommA_mean), stat="identity", width=0.01) +
  geom_errorbar(aes(xmax=CommA_mean+CommA_se, xmin=CommA_mean-CommA_se,y=CommB_mean), stat="identity", width=0.01) +
  labs(x = "CommA", y = "CommB") + geom_text(size=3) + theme_classic() + scale_color_manual(values=c("#393D47", "#990099"))


##CommBvsCommC
ggplot(data=wide, aes(x=CommC_mean, 
                      y=CommB_mean, color = CommC_day, label=CommC_variable))+ geom_point()+
  geom_smooth(method="lm",aes(x=CommC_mean,y=CommB_mean))+
  geom_errorbar(aes(ymax=CommB_mean+CommB_se, ymin=CommB_mean-CommB_se,x=CommC_mean), stat="identity", width=0.01) +
  geom_errorbar(aes(xmax=CommC_mean+CommC_se, xmin=CommC_mean-CommC_se,y=CommB_mean), stat="identity", width=0.01) +
  labs(x = "CommC", y = "CommB") + geom_text(size=3) + theme_classic()+ scale_color_manual(values=c("#393D47", "#990099"))

#### GLMMs ####
#relevel the predictor variables to ACM
data_filt <- data_filt %>% arrange(desc(treatment))
data_filt$treatment <- relevel(data_filt$treatment, ref = "ACM")

mbiomass <- brm(drymass_leaf1_grams ~ treatment,data=data_filt, family=Gamma(link="log"), iter = 10000, chains = 4, cores = 4)
saveRDS(mbiomass, file = "brms_mbiomass_ACMlevel.RDS")
mbiomass <- readRDS("brms_mbiomass_ACMlevel.RDS")
summary(mbiomass)
pp_check(mbiomass)


pmbiomass <- ggeffects::ggpredict(mbiomass, "treatment")
colnames(pmbiomass)[1] <- "treatment"
colnames(pmbiomass)[2] <- "drymass_leaf1_grams"

#Directional Effects
posterior_mbiomass <- as.data.frame(mbiomass)
CommA <- posterior_mbiomass %>% filter(b_treatmentCommA >0)
nrow(CommA)/nrow(posterior_mbiomass) #the probability of direction 0.8104
CommB <- posterior_mbiomass %>% filter(b_treatmentCommB >0)
nrow(CommB)/nrow(posterior_mbiomass) #the probability of direction .99725
CommC <- posterior_mbiomass %>% filter(b_treatmentCommC >0)
nrow(CommC)/nrow(posterior_mbiomass) #the probability of direction 0.8432
WATER <- posterior_mbiomass %>% filter(b_treatmentWATER <0)
nrow(WATER)/nrow(posterior_mbiomass) #the probability of direction 0.8125

posterior_mbiomass_melt <- posterior_mbiomass[,1:5]
posterior_mbiomass_melt <- reshape2::melt(posterior_mbiomass_melt)


level_order <- c("b_treatmentWATER", "b_treatmentCommC", "b_treatmentCommB", "b_treatmentCommA", "b_Intercept") 
ggplot(posterior_mbiomass_melt, aes(x = value, y =factor(variable, level = level_order),
  fill = variable)) +
  stat_halfeye(alpha=.6) +
  geom_vline(aes(xintercept=0), 
  color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + xlab("")+
  scale_fill_manual(values=c("#5C4033","#004488", "#ffaf49", "#44b7c2", "gray"))

#Fig3B remove intercept
posterior_mbiomass_melt_filt <- subset(posterior_mbiomass_melt, variable != "b_Intercept")
ggplot(posterior_mbiomass_melt_filt, aes(x = value, y =factor(variable, level = level_order),
                                    fill = variable)) +
  stat_halfeye(alpha=.6) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + xlab("")+
  scale_fill_manual(values=c("#004488", "#ffaf49", "#44b7c2", "gray"))+
  theme(legend.position = "none")



level_order2 <- c("CommA", "CommB", "CommC", "ACM", "WATER") 

#Fig2A
ggplot(data_filt, 
       mapping = aes(x = factor(treatment, level=level_order2), 
                     y = drymass_leaf1_grams, 
                     fill = treatment, color=treatment))  +
  geom_boxplot(outlier.shape=NA, alpha=.5) + geom_jitter(width=.4, size=3)+
  scale_color_manual(values=c("#5C4033","#004488", "#ffaf49", "#44b7c2", "gray")) +
  scale_fill_manual(values=c("#5C4033","#004488", "#ffaf49", "#44b7c2", "gray")) +
  geom_point(data=pmbiomass, aes(x=treatment, y=drymass_leaf1_grams), color="black", size=5) +
  geom_linerange(data=pmbiomass,aes(ymin=conf.low, ymax=conf.high), size=2, color="black",
                 position=position_dodge(width = 0.5)) + labs(
                   x = "Treatment", 
                   y = "Pitcher Dry Biomass (g)") + theme_classic()+
  theme(legend.position="none")+theme(
    axis.text.x = element_text(color = "black", size = 24),
    axis.text.y = element_text(color = "black", size = 24),
    axis.title.x = element_text(color = "black", size = 24),
    axis.title.y = element_text(color = "black", size = 24)
  )

#### leaf nutrients ~ treatment ####
m2 <- brm(nitrogen1_content_grams ~ treatment, data=data_filt, family=Gamma(link="log"), iter = 10000, chains = 4, cores = 4) 
m4 <- brm(carbon1_content_grams ~ treatment, data=data_filt, family=Gamma(link="log"), iter = 10000, chains = 4, cores = 4) 

saveRDS(m2, "leafN_treatment.RDS")
saveRDS(m4, "leafC_treatment.RDS")

m2 <- readRDS("leafN_treatment.RDS")
m4 <- readRDS("leafC_treatment.RDS")

summary(m2)
p2 <- ggpredict(m2, terms="treatment")

summary(m4)
p4 <- ggpredict(m4, terms="treatment")




#Directional Effect of N~treatment
posterior_m2 <- as.data.frame(m2)
CommA <- posterior_m2 %>% filter(b_treatmentCommA >0)
nrow(CommA)/nrow(posterior_m2)
CommB <- posterior_m2 %>% filter(b_treatmentCommB >0)
nrow(CommB)/nrow(posterior_m2)
CommC <- posterior_m2 %>% filter(b_treatmentCommC >0)
nrow(CommC)/nrow(posterior_m2)

WATER <- posterior_m2 %>% filter(b_treatmentWATER <0)
nrow(WATER)/nrow(posterior_m2)

posterior_interval(m2, prob=.95)
posterior_interval(m1, prob=.95)
posterior_interval(m3, prob=.95)

#marginal effect plots for C and N ~ treatment
pm2 <- ggpredict(m2, terms= "treatment")
colnames(pm2)[1] <- "treatment"
colnames(pm2)[2] <- "nitrogen1_content_grams"

#Fig2C
ggplot(data_filt, 
       mapping = aes(x = factor(treatment, level=level_order2), 
                     y = nitrogen1_content_grams, 
                     fill = treatment, color=treatment))  +
  geom_boxplot(outlier.shape=NA,alpha=.5) + geom_jitter(width=.4, size=3)+
  scale_color_manual(values=c("#5C4033","#004488", "#ffaf49", "#44b7c2", "gray")) +
  scale_fill_manual(values=c("#5C4033","#004488", "#ffaf49", "#44b7c2", "gray")) +
  geom_point(data=pm2, aes(x=treatment, y=nitrogen1_content_grams), color="black", size=5) +
  geom_linerange(data=pm2,aes(ymin=conf.low, ymax=conf.high), size=2, color="black",
                 position=position_dodge(width = 0.5)) + labs(
                   x = "Treatment", 
                   y = "Pitcher Nitrogen Content (g)") + theme_classic()+
  theme(legend.position="none")+theme(
    axis.text.x = element_text(color = "black", size = 24),
    axis.text.y = element_text(color = "black", size = 24),
    axis.title.x = element_text(color = "black", size = 24),
    axis.title.y = element_text(color = "black", size = 24)
  )

posterior_ncontent <- as.data.frame(m2)
posterior_ncontent_melt <- posterior_ncontent[,1:5]
posterior_ncontent_melt <- reshape2::melt(posterior_ncontent_melt)
level_order <- c("b_treatmentWATER", "b_treatmentCommC", "b_treatmentCommB", "b_treatmentCommA", "b_Intercept") 

#Fig2D
ggplot(posterior_ncontent_melt, aes(x = value, y =factor(variable, level = level_order),
                                    fill = variable)) +
  stat_halfeye(alpha=.6) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + xlab("")+
  scale_fill_manual(values=c("#5C4033","#004488", "#ffaf49", "#44b7c2", "gray"))+ theme(legend.position = "none")

#Fig 3D, intercept removed
posterior_ncontent_melt_filt <- subset(posterior_ncontent_melt, variable !="b_Intercept")

ggplot(posterior_ncontent_melt_filt, aes(x = value, y =factor(variable, level = level_order),
                                    fill = variable)) +
  stat_halfeye(alpha=.6) +
  geom_vline(aes(xintercept=0), 
             color="black", size=1, linetype="dashed")+
  ylab("Probability density") +
  theme_classic() + xlab("")+
  scale_fill_manual(values=c("#004488", "#ffaf49", "#44b7c2", "gray"))+ theme(legend.position = "none")




#m4, carbon~treatment
posteriorm4 <- mcmc_intervals_data(m4, 
                                   prob_outer=0.95,
                                   prob=0.5)

posteriorm4$nonzero <- NA
posteriorm4$nonzero[posteriorm4$ll>0 & posteriorm4$hh>0] <- "nonzero"
posteriorm4$nonzero[posteriorm4$ll<0 & posteriorm4$hh<0] <- "nonzero"
posteriorm4$nonzero[is.na(posteriorm4$nonzero)] <- "zero"
posteriorm4<- posteriorm4[1:5,]

level_order <- c("b_treatmentWATER","b_treatmentCommC", "b_treatmentCommB","b_treatmentCommA","b_Intercept") 
posteriorm4$parameter <- factor(posteriorm4$parameter, level = level_order)

#FIGS4D
ggplot(posteriorm4, aes(x = parameter,
                        shape=nonzero)) +
  geom_hline(yintercept = 0, linetype = 3, 
             size=1, color = "#b0b5b3") +
  geom_pointrange(aes(ymin = ll, ymax = hh, y = m),
                  position= position_dodge(width=0.75),
                  size = 3/4) +
  scale_color_manual(name="",
                     values = c("grey60", "#484c8d")) +
  scale_shape_manual(values=c(16, 17), 
                     labels=c("95% CI does\nnot contain zero", 
                              "95% CI\ncontains zero"))+
  coord_flip() +theme_bw() + 
  theme(axis.text.y = element_text( size=7), 
        axis.text.x=element_text(size=7),
        axis.title = element_text(size=7), 
        legend.text = element_text(size=7)) +
  xlab(NULL) +
  ylab("Estimated effect on pitcher carbon content")+
  guides(linetype=FALSE) + theme(legend.position = "none")

