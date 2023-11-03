library(tidyverse)
library(car)
library(lme4)
library(merTools)
library(mgcv)

#Import data
df <- read_csv("Data/pred.heat.tol.csv")%>%
  mutate(numDead = numTotal - numSurv)%>% #Add a number dead column
  mutate(across(c(species, trial, chamber, bath), as.factor))%>%
  mutate(species = fct_relevel(species, c("Buenoa", "Indica", "Irrorata", "Copto", "Pachy"))) #Relevel for graphing


#Convert data to long form
df.long<-df%>%
  uncount(numDead)%>%
  mutate(surv = 0)%>%
  dplyr::select(-numSurv)%>%
  rbind(mutate(uncount(dplyr::select(df, - numDead), numSurv), surv = 1))%>%
  arrange(trial, bath, species)
  

#Write a function to compute binomial GLM's for all species

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
#irrGlm <- run.glm("Irrorata") #Doesn't work because of complete separation of data. A different method will be used below
pachyGlm <- run.glm("Pachy")


##Make new data frame to add each species' data to for graphing
allHeatDf <- data.frame(fit = NULL, fit_link = NULL, se_link = NULL, fit_resp = NULL,
                          upr = NULL, lwr = NULL, species = NULL, temp = NULL)

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
#Back transform
Indicaheat<-mutate(ndata,
                   fit_resp=ilink(fit_link),
                   upr=ilink(fit_link+(se_link)),
                   lwr=ilink(fit_link-(se_link)), 
                   species = "Indica")

#Add to master df 
allHeatDf <- rbind(allHeatDf,Indicaheat)


##Buenoa

#Make new data set for predictions
Buen.new <- df.long%>%
  filter(species == "Buenoa")

ndata <- with(Buen.new, data_frame(temp = seq(min(temp), max(temp),
                                               length = 1000)))

ndata<-add_column(ndata,fit=predict(buenGlm, newdata=ndata, type="response"))

#Get inverse link function for prediction
ilink<-buenGlm$family$linkinv

###Add fit and se.fit on the link scale
ndata<-bind_cols(ndata, setNames(as_tibble(predict(
  buenGlm, ndata, se.fit = TRUE)[1:2]),
  c('fit_link','se_link')))

#Back transform
Buenheat<-mutate(ndata,
                 fit_resp=ilink(fit_link),
                 upr=ilink(fit_link+(se_link)),
                 lwr=ilink(fit_link-(se_link)), 
                 species = "Buenoa")

#Add to the master df
allHeatDf <- rbind (allHeatDf, Buenheat)


#Coptotomus

#Make new data set for predictions
Copto.new <- df.long%>%
  filter(species == "Copto")

ndata <- with(Copto.new, data_frame(temp = seq(min(temp), max(temp),
                                                length = 1000)))

ndata<-add_column(ndata,fit=predict(coptoGlm, newdata=ndata, type="response"))

#Get inverse link function for prediction
ilink<-coptoGlm$family$linkinv

###Add fit and se.fit on the link scale
ndata<-bind_cols(ndata, setNames(as_tibble(predict(
  coptoGlm, ndata, se.fit = TRUE)[1:2]),
  c('fit_link','se_link')))

#Back transform
Coptoheat<-mutate(ndata,
                  fit_resp=ilink(fit_link),
                  upr=ilink(fit_link+(se_link)),
                  lwr=ilink(fit_link-(se_link)), 
                  species = "Copto")

#Add to the master df
allHeatDf <- rbind (allHeatDf, Coptoheat)


#Pachy

#Make new data set for predictions
Pachy.new <- df.long%>%
  filter(species == "Pachy")

ndata <- with(Pachy.new, data_frame(temp = seq(min(temp), max(temp),
                                                length = 1000)))

ndata<-add_column(ndata,fit=predict(pachyGlm, newdata=ndata, type="response"))

#Get inverse link function for prediction
ilink<-pachyGlm$family$linkinv

###Add fit and se.fit on the link scale
ndata<-bind_cols(ndata, setNames(as_tibble(predict(
  pachyGlm, ndata, se.fit = TRUE)[1:2]),
  c('fit_link','se_link')))

#Back transform
Pachyheat<-mutate(ndata,
                  fit_resp=ilink(fit_link),
                  upr=ilink(fit_link+(se_link)),
                  lwr=ilink(fit_link-(se_link)), 
                  species = "Pachy")

#Add to the master df
allHeatDf <- rbind (allHeatDf, Pachyheat)


#Irrorata
#The Irrorata glm gave a warning because of complete/quasi separation of data. Conduct a Firth's biased-reduced logistic regression

#Subset data
Irrordat<-df.long%>%filter(species=="Irrorata")

#Conduct the Firth bias reduced log regression
library(logistf)
Irr.logistf<-logistf(formula=surv~temp, data=Irrordat, plcontrol = logistpl.control(maxit = 100000), control=logistf.control(maxstep=10, maxit=100000), pl=TRUE)

summary(Irr.logistf)

#Crate data set of predicted values
ndata <- with(Irrordat,data_frame(temp = seq(min(temp), max(temp),length = 1000)))
ndata<-add_column(ndata,fit = predict(Irr.logistf, newdata=ndata, type = "response"))




###Add fit and se.fit on the link scale
ndata<-bind_cols(ndata, setNames(as_tibble(predict(
  Irr.logistf, ndata, se.fit = TRUE)[1:2]),
  c('fit_link','se_link')))

#Back transform
Irrheat<-mutate(ndata,
                  fit_resp=ilink(fit_link),
                  upr=ilink(fit_link+(se_link)),
                  lwr=ilink(fit_link-(se_link)), 
                  species = "Irrorata")


#Add these predicted values to df with other species

allHeatDf <- rbind(allHeatDf, Irrheat)


#Calculate LT50 and add to a separate df
library(MASS)

#There is a function for generating the LT50's below, but it stopped working. Here is a less elegant method for each species individually

#Pachy
df <- data.frame(LT50 = NULL)  
for(i in 1:1000){
  
  n <- df.long%>%
    filter(species == "Pachy")%>%
    nrow()
  
  rand.df <- df.long%>%
    filter(species == "Pachy")%>% 
    sample_n(size = n, replace = T)%>%
    mutate(sim = i)
  
  model <- rand.df%>%
    filter(sim == i)%>%
    glm(data = ., surv ~ temp, family = "binomial")
  
  df <- data.frame(LT50 = unname(dose.p(model))[1],
                   species = "Pachy")%>%
    rbind(df)
  
}
PachyLT <- df

#Buenoa
df <- data.frame(LT50 = NULL)
for(i in 1:1000){
  
  n <- df.long%>%
    filter(species == "Buenoa")%>%
    nrow()
  
  rand.df <- df.long%>%
    filter(species == "Buenoa")%>% 
    sample_n(size = n, replace = T)%>%
    mutate(sim = i)
  
  model <- rand.df%>%
    filter(sim == i)%>%
    glm(data = ., surv ~ temp, family = "binomial")
  
  df <- data.frame(LT50 = unname(dose.p(model))[1],
                   species = "Buenoa")%>%
    rbind(df)
  
}
BuenLT <- df

#Copto
df <- data.frame(LT50 = NULL) 
for(i in 1:1000){
  
  n <- df.long%>%
    filter(species == "Copto")%>%
    nrow()
  
  rand.df <- df.long%>%
    filter(species == "Copto")%>% 
    sample_n(size = n, replace = T)%>%
    mutate(sim = i)
  
  model <- rand.df%>%
    filter(sim == i)%>%
    glm(data = ., surv ~ temp, family = "binomial")
  
  df <- data.frame(LT50 = unname(dose.p(model))[1],
                   species = "Copto")%>%
    rbind(df)
  
}
CoptLT <- df

#Indica
df <- data.frame(LT50 = NULL)
for(i in 1:1000){
  
  n <- df.long%>%
    filter(species == "Indica")%>%
    nrow()
  
  rand.df <- df.long%>%
    filter(species == "Indica")%>% 
    sample_n(size = n, replace = T)%>%
    mutate(sim = i)
  
  model <- rand.df%>%
    filter(sim == i)%>%
    glm(data = ., surv ~ temp, family = "binomial")
  
  df <- data.frame(LT50 = unname(dose.p(model))[1],
                   species = "Indica")%>%
    rbind(df)
  
}
IndLT <- df

#Bootstrap LT50 for Irrorata

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
  facet_wrap(~species, labeller = labeller(species = c("Buenoa" = "B. scimitra", 
                                                       "Indica" = "N. indica", 
                                                       "Irrorata" = "N. irrorata", 
                                                       "Copto" = "C. loticus", 
                                                       "Pachy" = "P. longipennis",
                                                       "Tramea" = "T. carolina"))) +
  geom_point(data = df.long, aes(y = surv), size = 2) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3, linetype = 0) +
  theme_classic() +
  labs(x = "Temperature (°C)", y = "Survivorship") +
  theme(axis.title = element_text(size = 28),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size = 24),
        strip.text = element_text(size=22, face = "italic"),
        legend.position = 0) +
  geom_vline(xintercept = 35, linetype = "dashed", size = 0.5, alpha = 0.7) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  geom_point(data = LT.avg.df, aes(x = temp, y = .5), color = "black", size = 2) + 
  geom_errorbarh(data = LT.avg.df, aes(xmin = lwrCI, xmax = uprCI, y = .5), color = "black", height = 0.05) 


#Uncomment to save
ggsave("Figures/Heat.tolerance.jpeg", width = 13.32, height = 7.27)














#####################################
#Here's the function method for generating LT50 that stopped working. It started underestimating the values for several species

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
