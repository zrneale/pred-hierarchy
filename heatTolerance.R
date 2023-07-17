
library(tidyverse)
library(car)
library(lme4)
library(merTools)
library(mgcv)

df <- read_csv("Data/pred.heat.tol.csv")%>%
  mutate(numDead = numTotal - numSurv)%>%
  mutate(across(c(species, trial, chamber, bath), as.factor))%>%
  mutate(species = fct_relevel(species, c("Buenoa", "Indica", "Irrorata", "Copto", "Pachy"))) #Relevel for graphing


#Convert to long form

df.long<-df%>%
  uncount(numDead)%>%
  mutate(surv = 0)%>%
  dplyr::select(-numSurv)%>%
  rbind(mutate(uncount(dplyr::select(df, - numDead), numSurv), surv = 1))%>%
  arrange(trial, bath, species)
  

#Write a function to compute binomial GLM's for all species

##Add columns for predicted values and SE limits to data frame
###Actually, that's going to take more time than I have. I'll do the less tidy way for now and tidy it up later

df.long%>%
  mutate(predict = NA, upperSE = NA, lowerSE = NA)

run.glm <- function(sp){
  #Calculate GLM
  model <- df.long%>%
  mutate(species = as.factor(species))%>%
  filter(species == sp)%>%
  glm(data = .,surv ~ temp, family = "binomial")
  

  return(model)
  
}

##Run the function for each species

buenGlm <- run.glm("Buenoa")
coptoGlm <- run.glm("Copto")
indGlm <- run.glm("Indica")
#irrGlm <- run.glm("Irrorata") #Doesn't work. I use a different method below
pachyGlm <- run.glm("Pachy")


##Make new data frame to add each species' data to for graphing
allHeatDf <- data.frame(fit = NA, fit_link = NA, se_link = NA, fit_resp = NA,
                          upr = NA, lwr = NA, species = NA, temp = NA)

#Calculate the predicted values with standard errors for each species


##make new data set for predictions
indDf <- df.long%>%
  filter(species == "Indica")

ndata <- with(indDf, data_frame(temp = seq(min(temp), max(temp),
                                                 length = 1000)))

ndata<-add_column(ndata,fit=predict(indGlm, newdata=ndata, type="response"))

##get inverse link function for prediction
ilink<-indGlm$family$linkinv

###add fit and se.fit on the link scale
ndata<-bind_cols(ndata, setNames(as_tibble(predict(
  indGlm, ndata, se.fit = TRUE)[1:2]),
  c('fit_link','se_link')))
#create tge ubtercak and back transform
Indicaheat<-mutate(ndata,
                   fit_resp=ilink(fit_link),
                   upr=ilink(fit_link+(se_link)),
                   lwr=ilink(fit_link-(se_link)), 
                   species = "Indica")

#get lethal temp 50 (50% die)
#lT50E<-ndata%>%mutate(fitR=round(fit, digits = 3))%>%filter(fitR>.4999 & fitR<0.501)
#lT50<-lT50E$Temp2



#plot Indicadata to check
Indicaheat%>%
  ggplot(aes(x = temp, y = fit)) +
  #geom_jitter(width = 0.5, height = 0) +
  geom_line(color="blue", size=3)+ 
  theme_classic() +
  labs(x = "Temperature (C)", y = "survivorship +/- SE") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size=16),
        strip.text = element_text(size=16))+
  geom_ribbon(data=Indicaheat, 
              alpha = 0.3, 
              aes(ymin = lwr, ymax = upr))+
  geom_point(data=indDf,aes(x=temp,y=surv), size=2, alpha=0.5)

#Looks good, add to master df 
allHeatDf <- rbind(allHeatDf,Indicaheat)


##Buenoa


#make new data set for predictions
Buen.new <- df.long%>%
  filter(species == "Buenoa")

ndata <- with(Buen.new, data_frame(temp = seq(min(temp), max(temp),
                                               length = 1000)))

ndata<-add_column(ndata,fit=predict(buenGlm, newdata=ndata, type="response"))

#get invers link function for prediction
ilink<-buenGlm$family$linkinv
###add fit and se.fit on the link scale
ndata<-bind_cols(ndata, setNames(as_tibble(predict(
  buenGlm, ndata, se.fit = TRUE)[1:2]),
  c('fit_link','se_link')))
#create tge ubtercak and back transform
Buenheat<-mutate(ndata,
                 fit_resp=ilink(fit_link),
                 upr=ilink(fit_link+(se_link)),
                 lwr=ilink(fit_link-(se_link)), 
                 species = "Buenoa")
#Plot it to check to see if any of these methods for calculating CI's put them on the correct scale
Buenheat%>%
  ggplot(aes(x = temp, y = fit)) +
  #geom_jitter(width = 0.5, height = 0) +
  geom_line(color="blue", size=3)+ 
  theme_classic() +
  labs(x = "Temperature (C)", y = "survivorship +/- SE") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size=16),
        strip.text = element_text(size=16))+
  geom_ribbon(data=Buenheat, 
              alpha = 0.3, 
              aes(ymin = lwr, ymax = upr))+
  geom_point(data=Buen.new,inherit.aes = F, aes(x=temp,y=surv), size=2, alpha=0.5)

#Add to the master df
allHeatDf <- rbind (allHeatDf, Buenheat)



#Coptotomus


#make new data set for predictions
Copto.new <- df.long%>%
  filter(species == "Copto")

ndata <- with(Copto.new, data_frame(temp = seq(min(temp), max(temp),
                                                length = 1000)))

ndata<-add_column(ndata,fit=predict(coptoGlm, newdata=ndata, type="response"))

#get invers link function for prediction
ilink<-coptoGlm$family$linkinv
###add fit and se.fit on the link scale
ndata<-bind_cols(ndata, setNames(as_tibble(predict(
  coptoGlm, ndata, se.fit = TRUE)[1:2]),
  c('fit_link','se_link')))
#create tge ubtercak and back transform
Coptoheat<-mutate(ndata,
                  fit_resp=ilink(fit_link),
                  upr=ilink(fit_link+(se_link)),
                  lwr=ilink(fit_link-(se_link)), 
                  species = "Copto")
#Plot it to check to see if any of these methods for calculating CI's put them on the correct scale
Coptoheat%>%
  ggplot(aes(x = temp, y = fit)) +
  #geom_jitter(width = 0.5, height = 0) +
  geom_line(color="blue", size=3)+ 
  theme_classic() +
  labs(x = "Temperature (C)", y = "survivorship +/- SE") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size=16),
        strip.text = element_text(size=16))+
  geom_ribbon(data=Coptoheat, 
              alpha = 0.3, 
              aes(ymin = lwr, ymax = upr))+
  geom_point(data=Copto.new,aes(x=temp,y=surv), size=2, alpha=0.5)

#Add to the master df
allHeatDf <- rbind (allHeatDf, Coptoheat)



#Pachy

#make new data set for predictions
Pachy.new <- df.long%>%
  filter(species == "Pachy")

ndata <- with(Pachy.new, data_frame(temp = seq(min(temp), max(temp),
                                                length = 1000)))

ndata<-add_column(ndata,fit=predict(pachyGlm, newdata=ndata, type="response"))

#get invers link function for prediction
ilink<-pachyGlm$family$linkinv
###add fit and se.fit on the link scale
ndata<-bind_cols(ndata, setNames(as_tibble(predict(
  pachyGlm, ndata, se.fit = TRUE)[1:2]),
  c('fit_link','se_link')))
#create tge ubtercak and back transform
Pachyheat<-mutate(ndata,
                  fit_resp=ilink(fit_link),
                  upr=ilink(fit_link+(se_link)),
                  lwr=ilink(fit_link-(se_link)), 
                  species = "Pachy")
#Plot it to check to see if any of these methods for calculating CI's put them on the correct scale
Pachyheat%>%
  ggplot(aes(x = temp, y = fit)) +
  #geom_jitter(width = 0.5, height = 0) +
  geom_line(color="blue", size=3)+ 
  theme_classic() +
  labs(x = "Temperature (C)", y = "survivorship +/- SE") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size=16),
        strip.text = element_text(size=16))+
  geom_ribbon(data=Pachyheat, 
              alpha = 0.3, 
              aes(ymin = lwr, ymax = upr))+
  geom_point(data=Pachy.new,aes(x=temp,y=surv), size=2, alpha=0.5)

#Add to the master df
allHeatDf <- rbind (allHeatDf, Pachyheat)


#The 4 easier ones are done. Now to do Irrorata. Volker sent some code that should work

##Irrorata

Irrordat<-df.long%>%filter(species=="Irrorata")
irrGlm<-glm(surv~temp, data=Irrordat, family="binomial",control=glm.control(maxit=500))
summary(irrGlm)
# issue is that we have complete /quasi complete separation, so there is only one temperature where we have 0 and 1 Temp2=97, everytihng <97 =1 and everything >97=0. So there is no variance to be explained
#Firth bias reduced log regression
library(logistf)
Irr.logistf<-logistf(formula=surv~temp, data=Irrordat, plcontrol = logistpl.control(maxit = 100000), control=logistf.control(maxstep=10, maxit=100000), pl=TRUE)
summary(Irr.logistf)

ndata <- with(Irrordat,data_frame(temp = seq(min(temp), max(temp),length = 100)))
ndata<-add_column(ndata,fit = predict(Irr.logistf, newdata=ndata, type = "response"))

ggplot(ndata, aes(x=temp, y=fit))+
  geom_vline(xintercept=35, linetype="dashed",size=2, color="grey")+
  geom_line(size=2, color="blue")+
  theme_classic()+
  geom_point(position=position_jitter(width=0.05,height=0.0), data=Irrordat, aes(x=temp, y=surv))

#That worked. Gonna need to add some columns to this data set in order to add it to the master df

allHeatDf <- ndata%>%
  mutate(fit_link = NA,
         se_link = NA,
         fit_resp = NA,
         upr = fit,
         lwr = fit,
         species = "Irrorata")%>%
  rbind(allHeatDf)


#Calculate LT50 and add to a separate df
library(MASS)


#Function to run the bootstrap 

getLT <- function(species, numsim){    
  
  df <- data.frame(LT50 = NULL)  
  
  for(i in 1:numsim){
    
    n <- df.long%>%
      filter(species == species)%>%
      nrow()
    
    rand.df <- df.long%>%
      filter(species == species)%>% 
      sample_n(size = n, replace = T)%>%
      mutate(sim = i)
    
    model <- rand.df%>%
      filter(sim == i)%>%
      glm(data = ., surv ~ temp, family = "binomial")
    
    df <- data.frame(LT50 = unname(dose.p(model))[1],
                     species = species)%>%
      rbind(df)
    
    
    
  }
  return(df)
}

#Run the function for Buenoa, Indica, Copto, and Pachy. I'll do a separate for loop for Irrorata since it uses a different model function

BuenLT <- getLT("Buenoa", 1000)
IndLT <- getLT("Indica", 1000)
CoptLT <- getLT("Copto", 1000)
PachyLT <- getLT("Pachy", 1000)

#Bootstrap for Irrorata

IrrLT <- data.frame(LT50 = NULL)

for(i in 1:1000){
  
  #Determine the sample size for the species
  n <- df.long%>%
    filter(species == "Irrorata")%>%
    nrow()
  
  #Create data set of resampled data
  rand.df <- df.long%>%
    filter(species == "Irrorata")%>% 
    sample_n(size = n, replace = T)%>%
    mutate(sim = i)
  
  #Run model with resampled data
  model <- rand.df%>%
    filter(sim == i)%>%
    logistf(formula = surv ~ temp, data = ., plcontrol = logistpl.control(maxit = 100000), control=logistf.control(maxstep=10, maxit=100000), pl=TRUE)
  
  IrrLT <- data.frame(LT50 = unname(dose.p(model))[1],
                   species = "Irrorata")%>%
    rbind(IrrLT)

}


#Create dataset with all species
LTdf <- data.frame(rbind(BuenLT,
                         IndLT,
                         IrrLT,
                         CoptLT,
                         PachyLT))%>%
  mutate(species = factor(species))



#Find average LT50 and confidence intervals
LT.avg.df <- LTdf%>%
  group_by(species)%>%
  dplyr::summarise(temp = mean(LT50),
                   lwrCI = quantile(LT50, 0.025),
                   uprCI = quantile(LT50, 0.975))



#Graph all species

##Change species in allHeatDf to factor and relevel for graphing
allHeatDf <- mutate_at(allHeatDf, vars(species), factor)%>%
  mutate(species = fct_relevel(species, c("Buenoa", "Indica", "Irrorata", "Copto", "Pachy")))

#Colorblind friendly palette
cbPalette <- c("#CC79A7", "#78C1EA", "#009E73", "#E69F00", "#D55E00", "#0072B2")

allHeatDf%>%
  drop_na(species)%>%
  ggplot(aes(x = temp, y = fit, color = species, fill = species)) +
  geom_line(size=1) +
  facet_wrap(~species, labeller = labeller(species = c("Buenoa" = "Buenoa sp.", 
                                                       "Indica" = "N. indica", 
                                                       "Irrorata" = "N. irrorata", 
                                                       "Copto" = "C. loticus", 
                                                       "Pachy" = "P. longipennis",
                                                       "Tramea" = "T. carolina"))) +
  geom_point(data = df.long, aes(y = surv), size = 2) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3, linetype = 0) +
  theme_classic() +
  labs(x = "Temperature (C)", y = "survivorship +/- SE") +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size=18),
        strip.text = element_text(size=16, face = "italic"),
        legend.position = 0) +
  #geom_vline(xintercept = LTdf$temp, linetype = "dotted", size = 2, alpha = 0.7) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  geom_point(data = LT.avg.df, aes(x = temp, y = .5), color = "black", size = 2) + 
  geom_errorbarh(data = LT.avg.df, aes(xmin = lwrCI, xmax = uprCI, y = .5), color = "black", height = 0.05) 


ggsave("Figures/Heat.tolerance.pdf", width = 13.32, height = 7.27)














#####################################


#Here's a couple other methods for calculating the LT50. I might want to go with one of these

##Trying a function that calculates LD50 using a logistic regression


library(HelpersMG)

#First I'll try with Buenoa
df%>%
  filter(species == "Buenoa")%>%
  data.frame()%>%
  mutate(alive = numSurv, dead = numDead, N = numTotal, doses = temp)%>%
  LD50()



#Trying a different function too, just to compare
library(ecotox)


df.long%>%
  filter(species == "Buenoa")%>%
  LC_logit(surv ~ temp, data = ., p = 50)


#Trying the dose.p function in Mass
##Write a function to create a data frame of all species

getLT <- function(species, model){
  data.frame(species = species,
             temp = unname(dose.p(model))[1],
             SE = dose.p(model)%>%attr("SE")%>%unname())
  
}

###Run function for each species and crate df.

LTdf <- rbind(getLT("Buenoa", buenGlm),
              getLT("Indica", indGlm),
              getLT("Irrorata", Irr.logistf),
              getLT("Copto", coptoGlm),
              getLT("Pachy", pachyGlm))%>%
  mutate(species = fct_relevel(species, c("Buenoa", "Indica", "Irrorata", "Copto", "Pachy"))) #Relevel for graphing


