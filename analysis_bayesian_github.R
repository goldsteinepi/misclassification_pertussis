#######################
# Analysis bayesian outcome correction
# Citation: Goldstein ND, Burstyn I, Newbern EC, Tabb LP, Gutowski J, Welles SL. Bayesian Correction of Misclassification of Pertussis in Vaccine Effectiveness Studies: How Much Does Underreporting Matter?
# 8/15/2014 -- Neal Goldstein
#######################


### FUNCTIONS ###

library(psych) #describe, describeBy
library(gmodels) #CrossTable
library(rjags)


### READ DATA ###

load("study_population_analysis.RData")
#load("study_population_analysis_CDRpotentials.RData")


### DROP N/A CONTROLS ###

matchid = studypop_final$matchid[1]
potentialcase = studypop_final$case[1]
case_asis = studypop_final$case[1]
studypop_final$control_sx = 0

#drop all controls from potential cases that were not cases (i.e., they are already controls)
for (i in 1:nrow(studypop_final))
{
  cat("\n\n**************** Observation: ",i," ****************\n",sep="")
  
  if ((potentialcase==1) && (case_asis==0))
  {
    if (matchid!=studypop_final$studyid[i]) {
      studypop_final$case_asis[i] = -1
      } else {
      studypop_final$control_sx[i] = 1
      }
  }
  
  if ((matchid != studypop_final$matchid[i+1]) && (!is.na(studypop_final$matchid[i+1])))
  {
    matchid = studypop_final$matchid[i+1]
    potentialcase = studypop_final$case[i+1]
    case_asis = studypop_final$case_asis[i+1]
  }
}
rm(i,matchid,case_asis,potentialcase)
studypop_final = studypop_final[studypop_final$case_asis!=-1 | is.na(studypop_final$case_asis),]
rownames(studypop_final) = NULL


### MERGE CDR POTENTIAL ###

#these potential (but not confirmed/probable) cases did not migrate from CDR to CDMS so had to manually join them
# studypop_final_CDRpotentials$control_sx = 1
# studypop_final = rbind(studypop_final, studypop_final_CDRpotentials)
# rm(studypop_final_CDRpotentials)


### SUBSET by CASE AGE for UTD ###

population = NA
matchid = studypop_final$matchid[1]
case_age = studypop_final$age_weeks[1]

#only include cases that can be classified as UTD/NUTD (min age 13 weeks)
for (i in 1:nrow(studypop_final))
{
  cat("\n\n**************** Observation: ",i," ****************\n",sep="")
  
  if (case_age>=13)
  #if ((case_age>=13) && (case_age<52)) #limit analysis to infants only <52 weeks (this is specific to cases, controls may be up to 2 weeks older)
    population = rbind(population, studypop_final[i, ])
  
  if ((matchid != studypop_final$matchid[i+1]) && (!is.na(studypop_final$matchid[i+1])))
  {
    matchid = studypop_final$matchid[i+1]
    case_age = studypop_final$age_weeks[i+1]
  }
}
rm(i,matchid,case_age)
population = population[-1,]
rownames(population) = NULL


### RECODE ###

#collapse probable and confirmed
population$case_asis = ifelse(population$case_asis==0 | is.na(population$case_asis), 0, 1)
population$case_1997 = ifelse(population$case_1997==0 | is.na(population$case_1997), 0, 1)
population$case_2014 = ifelse(population$case_2014==0 | is.na(population$case_2014), 0, 1)

#dichotomize parity
population$mother_parous = ifelse(population$mother_parity>=1, 1, 0)

#dichotomize age
population$age_1yr = ifelse(population$age==0,0,1)

#dichtomize race
population$race_black = ifelse(population$race==1,1,0)

#doses indicator variables
population$dose0 = ifelse(population$pertussis_vax==0, 1, 0)
population$dose1 = ifelse(population$pertussis_vax>=1, 1, 0)
population$dose2 = ifelse(population$pertussis_vax>=2, 1, 0)
population$dose3 = ifelse(population$pertussis_vax>=3, 1, 0)
population$dose4 = ifelse(population$pertussis_vax>=4, 1, 0)
population$dose5 = ifelse(population$pertussis_vax>=5, 1, 0)


### SUBSET DATA ###

#cdms-era only
population = population[population$case_year>=2011, ]

# validation = NA
# analysis = NA
# dataset = "A"
# matchid = population$matchid[1]
# population$dataset = NA
# 
# #create two subsets of data: 1 validation & 1 analysis
# for (i in 1:nrow(population))
# {
#   cat("\n\n**************** Observation: ",i," ****************\n",sep="")
# 
#   if (dataset=="V") {
#     validation = rbind(validation, population[i, ])
#     population$dataset[i] = "V"
#   } else if (dataset=="A") {
#     analysis = rbind(analysis, population[i, ])
#     population$dataset[i] = "A"
#   }
#   
#   
#   if ((matchid != population$matchid[i+1]) && (!is.na(population$matchid[i+1])))
#   {
#     matchid = population$matchid[i+1]
#     
#     #flip flop dataset
#     if (dataset=="V") {
#       dataset = "A"
#     } else {
#       dataset = "V"
#     }
#   }
# }
# rm(i,dataset,matchid)
# analysis = analysis[-1,]
# validation = validation[-1,]
# rownames(analysis) = NULL
# rownames(validation) = NULL


### ANALYSIS: POPULATION CHARACTERISTICS ###

CrossTable(population$age_1yr)
CrossTable(population$gender)
CrossTable(population$race_black)
CrossTable(population$ethnicity)
describe(population$mother_age)
CrossTable(population$mother_married)
CrossTable(population$mother_foreign_born)
CrossTable(population$mother_education)
CrossTable(population$mother_insurance)
CrossTable(population$mother_parous)
CrossTable(population$pertussis_vax_utd)

#pertussis classification
CrossTable(population$case)
CrossTable(population$case_asis[population$case==1])


### ANALYSIS: COMPARISON of APPARENT CASE and HEALTH CONTROLS ###

CrossTable(population$age_1yr[population$control_sx==0], population$case_asis[population$control_sx==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$gender[population$control_sx==0], population$case_asis[population$control_sx==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$race_black[population$control_sx==0], population$case_asis[population$control_sx==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$ethnicity[population$control_sx==0], population$case_asis[population$control_sx==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(population$mother_age[population$control_sx==0], population$case_asis[population$control_sx==0]); t.test(population$mother_age[population$control_sx==0] ~ population$case_asis[population$control_sx==0])
CrossTable(population$mother_married[population$control_sx==0], population$case_asis[population$control_sx==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$mother_foreign_born[population$control_sx==0], population$case_asis[population$control_sx==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$mother_education[population$control_sx==0], population$case_asis[population$control_sx==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$mother_insurance[population$control_sx==0], population$case_asis[population$control_sx==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$mother_parous[population$control_sx==0], population$case_asis[population$control_sx==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$pertussis_vax_utd[population$control_sx==0], population$case_asis[population$control_sx==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)


### ANALYSIS: COMPARISON of Sx CONTROLS and HEALTHY CONTROLS ###

CrossTable(population$age_1yr[population$case_asis==0], population$control_sx[population$case_asis==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$gender[population$case_asis==0], population$control_sx[population$case_asis==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$race_black[population$case_asis==0], population$control_sx[population$case_asis==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$ethnicity[population$case_asis==0], population$control_sx[population$case_asis==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(population$mother_age[population$case_asis==0], population$control_sx[population$case_asis==0]); t.test(population$mother_age[population$case_asis==0] ~ population$control_sx[population$case_asis==0])
CrossTable(population$mother_married[population$case_asis==0], population$control_sx[population$case_asis==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$mother_foreign_born[population$case_asis==0], population$control_sx[population$case_asis==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$mother_education[population$case_asis==0], population$control_sx[population$case_asis==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$mother_insurance[population$case_asis==0], population$control_sx[population$case_asis==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$mother_parous[population$case_asis==0], population$control_sx[population$case_asis==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$pertussis_vax_utd[population$case_asis==0], population$control_sx[population$case_asis==0], prop.r=F, prop.t=F, prop.chisq=F, chisq=T)


# ### ANALYSIS: COMPARISON of 2011 VALIDATION and REST OF POPULATION ###
# 
# CrossTable(population$age_1yr, population$case_year>=2011, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$gender, population$case_year>=2011, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$race_black, population$case_year>=2011, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$ethnicity, population$case_year>=2011, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# describeBy(population$mother_age, population$case_year>=2011); t.test(population$mother_age ~ population$case_year>=2011)
# CrossTable(population$mother_married, population$case_year>=2011, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$mother_foreign_born, population$case_year>=2011, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$mother_education, population$case_year>=2011, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$mother_insurance, population$case_year>=2011, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$mother_parous, population$case_year>=2011, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$pertussis_vax_utd, population$case_year>=2011, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)


### ANALYSIS: Comparison by missingness ###

#missing marker
population$missing = ifelse(is.na(population$age_1yr) | is.na(population$race_black) | is.na(population$mother_parous) | is.na(population$pertussis_vax_utd), 1, 0)
CrossTable(population$missing)

CrossTable(population$age_1yr, population$missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$gender, population$missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$race_black, population$missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$ethnicity, population$missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(population$mother_age, population$missing); t.test(population$mother_age ~ population$missing)
CrossTable(population$mother_married, population$missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$mother_foreign_born, population$missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$mother_education, population$missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$mother_insurance, population$missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$mother_parous, population$missing, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)


### ANALYSIS: NAIVE ###

model = glm(case_asis~as.factor(pertussis_vax_utd)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=population,family=binomial(link="logit"))
summary(model)
round(exp(coef(model)),2)
round(exp(confint(model)),2)

#2x2 using misclassified and perfect
CrossTable(population$pertussis_vax_utd, population$case_asis, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
#CrossTable(population$pertussis_vax_utd, population$case_1997, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$pertussis_vax_utd, population$case_2014, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)


### PRIOR SPECIFICATIONS ###

# #nondifferential, investigator
# CrossTable(population$case_asis[population$case==1], population$case_1997[population$case==1], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$case_asis[population$case==1], population$case_2014[population$case==1], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# 
# #differential, investigator
# CrossTable(population$case_asis[population$case==1 & population$pertussis_vax_utd==0], population$case_1997[population$case==1 & population$pertussis_vax_utd==0], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$case_asis[population$case==1 & population$pertussis_vax_utd==0], population$case_2014[population$case==1 & population$pertussis_vax_utd==0], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$case_asis[population$case==1 & population$pertussis_vax_utd==1], population$case_1997[population$case==1 & population$pertussis_vax_utd==1], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$case_asis[population$case==1 & population$pertussis_vax_utd==1], population$case_2014[population$case==1 & population$pertussis_vax_utd==1], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)

#nondifferential, 2011 CDMS investigator
#CrossTable(population$case_asis[population$case_year>=2011 & population$case==1], population$case_1997[population$case_year>=2011 & population$case==1], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(population$case_asis[population$case_year>=2011 & population$case==1], population$case_2014[population$case_year>=2011 & population$case==1], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)

# #differential, 2011 CDMS investigator
# CrossTable(population$case_asis[population$case_year>=2011 & population$case==1 & population$pertussis_vax_utd==0], population$case_1997[population$case_year>=2011 & population$case==1 & population$pertussis_vax_utd==0], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$case_asis[population$case_year>=2011 & population$case==1 & population$pertussis_vax_utd==0], population$case_2014[population$case_year>=2011 & population$case==1 & population$pertussis_vax_utd==0], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$case_asis[population$case_year>=2011 & population$case==1 & population$pertussis_vax_utd==1], population$case_1997[population$case_year>=2011 & population$case==1 & population$pertussis_vax_utd==1], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$case_asis[population$case_year>=2011 & population$case==1 & population$pertussis_vax_utd==1], population$case_2014[population$case_year>=2011 & population$case==1 & population$pertussis_vax_utd==1], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# 
# #nondifferential, population
# CrossTable(population$case, population$case_1997, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$case, population$case_2014, prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# 
# #differential, population
# CrossTable(population$case[population$pertussis_vax_utd==0], population$case_1997[population$pertussis_vax_utd==0], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$case[population$pertussis_vax_utd==0], population$case_2014[population$pertussis_vax_utd==0], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$case[population$pertussis_vax_utd==1], population$case_1997[population$pertussis_vax_utd==1], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$case[population$pertussis_vax_utd==1], population$case_2014[population$pertussis_vax_utd==1], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# 
# #nondifferential, 2011 CDMS population
# CrossTable(population$case[population$case_year>=2011], population$case_1997[population$case_year>=2011], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$case[population$case_year>=2011], population$case_2014[population$case_year>=2011], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# 
# #differential, 2011 CDMS population
# CrossTable(population$case[population$case_year>=2011 & population$pertussis_vax_utd==0], population$case_1997[population$case_year>=2011 & population$pertussis_vax_utd==0], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$case[population$case_year>=2011 & population$pertussis_vax_utd==0], population$case_2014[population$case_year>=2011 & population$pertussis_vax_utd==0], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$case[population$case_year>=2011 & population$pertussis_vax_utd==1], population$case_1997[population$case_year>=2011 & population$pertussis_vax_utd==1], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(population$case[population$case_year>=2011 & population$pertussis_vax_utd==1], population$case_2014[population$case_year>=2011 & population$pertussis_vax_utd==1], prop.c=F, prop.t=F, prop.chisq=F, chisq=T)
# 
# #prior OR 1997, 2011 CDMS
# model = glm(case_1997~as.factor(pertussis_vax_utd)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=subset(population,population$case_year>=2011),family=binomial(link="logit"))
# summary(model)
# round(coef(model),2)
# #variance is standard error (apprxs the standard deviation) squared
# round(0.2585*0.2585,2)
# 
# round(exp(coef(model)),2)
# round(exp(confint(model)),2)
# 
# #prior OR 2014, 2011 CDMS
# model = glm(case_2014~as.factor(pertussis_vax_utd)+as.factor(age_1yr)+as.factor(race)+as.factor(mother_parous),data=subset(population,population$case_year>=2011),family=binomial(link="logit"))
# summary(model)
# round(coef(model),2)
# #variance is standard error (apprxs the standard deviation) squared
# round(0.25060*0.25060,2)
# 
# round(exp(coef(model)),2)
# round(exp(confint(model)),2)


### BUGS MODEL ###

#2011 nondifferential final
bugs_model =
  "model {
    for (i in 1:n) {
    
      #outcome model, log odds of pertussis given these predictors
      case[i] ~ dbern(p_case[i])
      logit(p_case[i]) <- b0+b1*utd[i]+b2*age[i]+b3*race[i]+b4*multiparous[i]
      
      #exposure models, log odds of utd given these predictors
      utd[i] ~ dbern(p_utd[i])
      logit(p_utd[i]) <- a0+a2*age[i]+a3*race[i]+a4*multiparous[i]
      
      #measurement model, imputing the true case status given the measurement error
      case_star[i] ~ dbern(p_case_star[i])
      p_case_star[i] <- sn0*case[i]+(1-case[i])*(1-sp0)
      
      #prevalence models of potential confounders
      age[i] ~ dbern(p_age[i])
      logit(p_age[i]) <- prev_age
      
      race[i] ~ dbern(p_race[i])
      logit(p_race[i]) <- prev_race
      
      multiparous[i] ~ dbern(p_multiparous[i])
      logit(p_multiparous[i]) <- prev_multiparous
    
    }
    
    #priors
    #for normal distribution, provide (mean, precision=(1/variance))
    #for beta distribution, provide (alpha, beta), add beta(1,1)
    
    b0 ~ dnorm(0,1/10)
    #b1 ~ dnorm(-1.77,1/0.002)
    b1 ~ dnorm(-1.77,1/0.5)
    b2 ~ dnorm(0,1/10)
    b3 ~ dnorm(0,1/10)
    b4 ~ dnorm(0,1/10)
    a0 ~ dnorm(0,1/10)
    a2 ~ dnorm(0,1/10)
    a3 ~ dnorm(0,1/10)
    a4 ~ dnorm(0,1/10)
    prev_age ~ dnorm(0,1/10)
    prev_race ~ dnorm(0,1/10)
    prev_multiparous ~ dnorm(0,1/10)
    sn0 ~ dbeta(113,13) #E=0
    sp0 ~ dbeta(119873,1) #E=0
    #sn0 ~ dbeta(113,18) #E=5
    #sp0 ~ dbeta(119868,1) #E=5
    #sn0 ~ dbeta(113,63) #E=50
    #sp0 ~ dbeta(119823,1) #E=50
    #sn0 ~ dbeta(113,513) #E=500
    #sp0 ~ dbeta(119373,1) #E=500
}"

#2011 CDMS investigator
bugs_model =
"model {
  for (i in 1:n) {
  
    #outcome model, log odds of pertussis given these predictors
    case[i] ~ dbern(p_case[i])
    logit(p_case[i]) <- b0+b1*utd[i]+b2*age[i]+b3*race[i]+b4*multiparous[i]
    
    #exposure models, log odds of utd given these predictors
    utd[i] ~ dbern(p_utd[i])
    logit(p_utd[i]) <- a0+a2*age[i]+a3*race[i]+a4*multiparous[i]
    
    #measurement model, imputing the true case status given the measurement error
    #this only applies to potential cases, otherwise they are a healthy control (SN=SP=1)
    case_star[i] ~ dbern(p_case_star[i])
    p_case_star[i] <- potential[i]*(sn0*case[i]*(1-utd[i])+(1-case[i])*(1-sp0)*(1-utd[i]) + sn1*case[i]*(utd[i])+(1-case[i])*(1-sp1)*(utd[i]))
    
    #prevalence models of potential confounders
    age[i] ~ dbern(p_age[i])
    logit(p_age[i]) <- prev_age
    
    race[i] ~ dbern(p_race[i])
    logit(p_race[i]) <- prev_race
    
    multiparous[i] ~ dbern(p_multiparous[i])
    logit(p_multiparous[i]) <- prev_multiparous
    
  }

  #priors
  #for normal distribution, provide (mean, precision=(1/variance))
  #for beta distribution, provide (alpha, beta), add beta(1,1)
  
  b0 ~ dnorm(0,1/10)
  #b1 ~ dnorm(-0.57,1/0.07) #1997
  b1 ~ dnorm(-0.65,1/0.06) #2014
  b2 ~ dnorm(0,1/10)
  b3 ~ dnorm(0,1/10)
  b4 ~ dnorm(0,1/10)
  a0 ~ dnorm(0,1/10)
  a2 ~ dnorm(0,1/10)
  a3 ~ dnorm(0,1/10)
  a4 ~ dnorm(0,1/10)
  prev_age ~ dnorm(0,1/10)
  prev_race ~ dnorm(0,1/10)
  prev_multiparous ~ dnorm(0,1/10)
  #sn0 ~ dbeta(39,1) #1997
  #sp0 ~ dbeta(51,2) #1997
  #sn1 ~ dbeta(75,8) #1997
  #sp1 ~ dbeta(179,2) #1997
  sn0 ~ dbeta(40,3) #2014
  sp0 ~ dbeta(49,1) #2014
  sn1 ~ dbeta(75,12) #2014
  sp1 ~ dbeta(175,2) #2014
}"

#2011 CDMS population
bugs_model =
"model {
  for (i in 1:n) {
  
    #outcome model, log odds of pertussis given these predictors
    case[i] ~ dbern(p_case[i])
    logit(p_case[i]) <- b0+b1*utd[i]+b2*age[i]+b3*race[i]+b4*multiparous[i]
    
    #exposure models, log odds of utd given these predictors
    utd[i] ~ dbern(p_utd[i])
    logit(p_utd[i]) <- a0+a2*age[i]+a3*race[i]+a4*multiparous[i]
    
    #measurement model, imputing the true case status given the measurement error
    #this only applies to potential cases, otherwise they are a healthy control (SN=SP=1)
    case_star[i] ~ dbern(p_case_star[i])
    p_case_star[i] <- sn0*case[i]*(1-utd[i])+(1-case[i])*(1-sp0)*(1-utd[i]) + sn1*case[i]*(utd[i])+(1-case[i])*(1-sp1)*(utd[i])
    
    #prevalence models of potential confounders
    age[i] ~ dbern(p_age[i])
    logit(p_age[i]) <- prev_age
    
    race[i] ~ dbern(p_race[i])
    logit(p_race[i]) <- prev_race
    
    multiparous[i] ~ dbern(p_multiparous[i])
    logit(p_multiparous[i]) <- prev_multiparous
  
  }
  
  #priors
  #for normal distribution, provide (mean, precision=(1/variance))
  #for beta distribution, provide (alpha, beta), add beta(1,1)
  
  b0 ~ dnorm(0,1/10)
  #b1 ~ dnorm(-0.57,1/0.07) #1997
  b1 ~ dnorm(-0.65,1/0.06) #2014
  b2 ~ dnorm(0,1/10)
  b3 ~ dnorm(0,1/10)
  b4 ~ dnorm(0,1/10)
  a0 ~ dnorm(0,1/10)
  a2 ~ dnorm(0,1/10)
  a3 ~ dnorm(0,1/10)
  a4 ~ dnorm(0,1/10)
  prev_age ~ dnorm(0,1/10)
  prev_race ~ dnorm(0,1/10)
  prev_multiparous ~ dnorm(0,1/10)
  #sn0 ~ dbeta(39,1) #1997
  #sp0 ~ dbeta(20,52) #1997
  #sn1 ~ dbeta(82,1) #1997
  #sp1 ~ dbeta(93,180) #1997
  sn0 ~ dbeta(42,1) #2014
  sp0 ~ dbeta(20,49) #2014
  sn1 ~ dbeta(86,1) #2014
  sp1 ~ dbeta(93,176) #2014
}"


### COMPLETE CASE ###

analysis_complete = na.omit(population[,c("case_asis","pertussis_vax_utd","age_1yr","race_black","mother_parous","case")])


### BAYESIAN SAMPLING ###

#initialize model
jags <- jags.model(textConnection(bugs_model),
                   data = list('case_star' = analysis_complete$case_asis,
                               'utd' = analysis_complete$pertussis_vax_utd,
                               'age' = analysis_complete$age_1yr,
                               'race' = analysis_complete$race_black,
                               'multiparous' = analysis_complete$mother_parous,
                               'n' = nrow(analysis_complete)),
                   n.chains = 2,
                   n.adapt = 100)

#sample from the posterior distribution
#jags_samples = coda.samples(jags, variable.names=c("b1","sn0","sp0","sn1","sp1","case"), n.iter=20000)
jags_samples = coda.samples(jags, variable.names=c("b1","sn0","sp0","case"), n.iter=20000)


### SAVE POSTERIORS ###

#E=0
#save.image(Posterior 0.RData")
#load(Posterior 0.RData")

#E=5
#save.image(Posterior 5.RData")
#load(Posterior 5.RData")

#E=50
#save.image(Posterior 50.RData")
#load(Posterior 50.RData")

#E=500
#save.image(Posterior 500.RData")
#load(Posterior 500.RData")

#1997 CDMS 2011 investigator
#save.image(1997 investigator.RData")
#load(1997 investigator.RData")

#2014 CDMS 2011 investigator
#save.image(2014 investigator.RData")
#load(2014 investigator.RData")

#1997 CDMS 2011 population
#save.image(1997 population.RData")
#load(1997 population.RData")

#2014 CDMS 2011 population
#save.image(2014 population.RData")
#load(2014 population.RData")


### BAYESIAN INFERENCE ###

#check for convergence
plot(jags_samples[,"b1"]) #good mixing and shape of distribution
gelman.plot(jags_samples[,"b1"]) #no discernable difference
#gelman.diag(jags_samples[,c("b1","sn0","sp0","sn1","sp1")]) #<1.1
gelman.diag(jags_samples[,c("b1","sn0","sp0")]) #<1.1

#statistics for model fit, discard first 1000 observations for burn in
#summary(window(jags_samples[,c("b1","sn0","sp0","sn1","sp1")], start=1000))
summary(window(jags_samples[,c("b1","sn0","sp0")], start=1000))


### COUNT NUMBER OF CASES ###

#can't tabulate a prevalance since this is a case control study

#obtain outcome status for each individual, for each simulation, in each chain
prev_chain1 = as.data.frame(window(jags_samples[[1]][,2:(nrow(analysis_complete)+1)], start=1000))
prev_chain2 = as.data.frame(window(jags_samples[[2]][,2:(nrow(analysis_complete)+1)], start=1000))

#tally the total number for prevalance
n_cases = c(rowSums(prev_chain1), rowSums(prev_chain2))

#summary of prevalence
mean(n_cases)
quantile(n_cases,probs=c(0.025,0.975))

