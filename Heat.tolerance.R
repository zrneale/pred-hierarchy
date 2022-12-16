
library(tidyverse)
library(car)
library(lme4)
library(merTools)
library(mgcv)

df <- read_csv("Data/Pred\ heat\ tolerance.csv")%>%
  mutate(Numdead = Totalnum - Numsurv)%>%
  mutate(Species = as.factor(Species))

#Convert to long form

df.long<-df%>%
  uncount(Numdead)%>%
  mutate(Surv = 0)%>%
  dplyr::select(-Numsurv)%>%
  rbind(mutate(uncount(dplyr::select(df, - Numdead), Numsurv), Surv = 1))%>%
  arrange(Trial, Bath, Species) %>%
  mutate(Species = fct_relevel(Species, c("Buenoa", "Indica", "Irrorata", "Copto", "Pachy"))) #Relevel for graphing


#Write a function to compute binomial GLM's for all species

##Add columns for predicteds and SE limits to data frame
###Actually, that's going to take more time than I have. I'll do the less tidy way for now and tidy it up later

df.long%>%
  mutate(Predict = NA, UpperSE = NA, LowerSE = NA)
run.glm <- function(sp){
  #Calculate GLM
  model <- df.long%>%
  mutate(Species = as.factor(Species))%>%
  filter(Species == sp)%>%
  glm(data = .,Surv ~ TempC, family = "binomial")
  

  return(model)
  
}

##Run the function for each species

Buen.glm <- run.glm("Buenoa")
Copto.glm <- run.glm("Copto")
Indica.glm <- run.glm("Indica")
#Irrorata.glm <- run.gam("Irrorata") Doesn't work. I use a different method below
Pachy.glm <- run.glm("Pachy")


#Try getting CI values that fit with 0,1 bounds with Indica alone. I'll also make a new data frame for all the species to be added together for the final graph

##Make new data frame to add each species' data to for graphing
Heat.all.df <- data.frame(fit = NA, fit_link = NA, se_link = NA, fit_resp = NA,
                          upr = NA, lwr = NA, Species = NA, TempC = NA)


##make new data set for predictions
Indica.new <- df.long%>%
  filter(Species == "Indica")

ndata <- with(Indica.new, data_frame(TempC = seq(min(TempC), max(TempC),
                                                 length = 1000)))

ndata<-add_column(ndata,fit=predict(Indica.glm, newdata=ndata, type="response"))

##get inverse link function for prediction
ilink<-Indica.glm$family$linkinv

###add fit and se.fit on the link scale
ndata<-bind_cols(ndata, setNames(as_tibble(predict(
  Indica.glm, ndata, se.fit = TRUE)[1:2]),
  c('fit_link','se_link')))
#create tge ubtercak and back transform
Indicaheat<-mutate(ndata,
                   fit_resp=ilink(fit_link),
                   upr=ilink(fit_link+(se_link)),
                   lwr=ilink(fit_link-(se_link)), 
                   Species = "Indica")

#get lethal temp 50 (50% die)
#lT50E<-ndata%>%mutate(fitR=round(fit, digits = 3))%>%filter(fitR>.4999 & fitR<0.501)
#lT50<-lT50E$Temp2



#plot Indicadata to check
Indicaheat%>%
  ggplot(aes(x = TempC, y = fit)) +
  #geom_jitter(width = 0.5, height = 0) +
  geom_line(color="blue", size=3)+ 
  theme_classic() +
  labs(x = "Temperature (C)", y = "Survivorship +/- SE") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size=16),
        strip.text = element_text(size=16))+
  geom_ribbon(data=Indicaheat, 
              alpha = 0.3, 
              aes(ymin = lwr, ymax = upr))+
  geom_point(data=Indica.new,aes(x=TempC,y=Surv), size=2, alpha=0.5)

#Looks good, add to master df 
Heat.all.df <- rbind(Heat.all.df,Indicaheat)


##Buenoa


#make new data set for predictions
Buen.new <- df.long%>%
  filter(Species == "Buenoa")

ndata <- with(Buen.new, data_frame(TempC = seq(min(TempC), max(TempC),
                                               length = 1000)))

ndata<-add_column(ndata,fit=predict(Buen.glm, newdata=ndata, type="response"))

#get invers link function for prediction
ilink<-Buen.glm$family$linkinv
###add fit and se.fit on the link scale
ndata<-bind_cols(ndata, setNames(as_tibble(predict(
  Buen.glm, ndata, se.fit = TRUE)[1:2]),
  c('fit_link','se_link')))
#create tge ubtercak and back transform
Buenheat<-mutate(ndata,
                 fit_resp=ilink(fit_link),
                 upr=ilink(fit_link+(se_link)),
                 lwr=ilink(fit_link-(se_link)), 
                 Species = "Buenoa")
#Plot it to check to see if any of these methods for calculating CI's put them on the correct scale
Buenheat%>%
  ggplot(aes(x = TempC, y = fit)) +
  #geom_jitter(width = 0.5, height = 0) +
  geom_line(color="blue", size=3)+ 
  theme_classic() +
  labs(x = "Temperature (C)", y = "Survivorship +/- SE") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size=16),
        strip.text = element_text(size=16))+
  geom_ribbon(data=Buenheat, 
              alpha = 0.3, 
              aes(ymin = lwr, ymax = upr))+
  geom_point(data=Buen.new,inherit.aes = F, aes(x=TempC,y=Surv), size=2, alpha=0.5)

#Add to the master df
Heat.all.df <- rbind (Heat.all.df, Buenheat)



#Coptotomus


#make new data set for predictions
Copto.new <- df.long%>%
  filter(Species == "Copto")

ndata <- with(Copto.new, data_frame(TempC = seq(min(TempC), max(TempC),
                                                length = 1000)))

ndata<-add_column(ndata,fit=predict(Copto.glm, newdata=ndata, type="response"))

#get invers link function for prediction
ilink<-Copto.glm$family$linkinv
###add fit and se.fit on the link scale
ndata<-bind_cols(ndata, setNames(as_tibble(predict(
  Copto.glm, ndata, se.fit = TRUE)[1:2]),
  c('fit_link','se_link')))
#create tge ubtercak and back transform
Coptoheat<-mutate(ndata,
                  fit_resp=ilink(fit_link),
                  upr=ilink(fit_link+(se_link)),
                  lwr=ilink(fit_link-(se_link)), 
                  Species = "Copto")
#Plot it to check to see if any of these methods for calculating CI's put them on the correct scale
Coptoheat%>%
  ggplot(aes(x = TempC, y = fit)) +
  #geom_jitter(width = 0.5, height = 0) +
  geom_line(color="blue", size=3)+ 
  theme_classic() +
  labs(x = "Temperature (C)", y = "Survivorship +/- SE") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size=16),
        strip.text = element_text(size=16))+
  geom_ribbon(data=Coptoheat, 
              alpha = 0.3, 
              aes(ymin = lwr, ymax = upr))+
  geom_point(data=Copto.new,aes(x=TempC,y=Surv), size=2, alpha=0.5)

#Add to the master df
Heat.all.df <- rbind (Heat.all.df, Coptoheat)



#Pachy

#make new data set for predictions
Pachy.new <- df.long%>%
  filter(Species == "Pachy")

ndata <- with(Pachy.new, data_frame(TempC = seq(min(TempC), max(TempC),
                                                length = 1000)))

ndata<-add_column(ndata,fit=predict(Pachy.glm, newdata=ndata, type="response"))

#get invers link function for prediction
ilink<-Pachy.glm$family$linkinv
###add fit and se.fit on the link scale
ndata<-bind_cols(ndata, setNames(as_tibble(predict(
  Pachy.glm, ndata, se.fit = TRUE)[1:2]),
  c('fit_link','se_link')))
#create tge ubtercak and back transform
Pachyheat<-mutate(ndata,
                  fit_resp=ilink(fit_link),
                  upr=ilink(fit_link+(se_link)),
                  lwr=ilink(fit_link-(se_link)), 
                  Species = "Pachy")
#Plot it to check to see if any of these methods for calculating CI's put them on the correct scale
Pachyheat%>%
  ggplot(aes(x = TempC, y = fit)) +
  #geom_jitter(width = 0.5, height = 0) +
  geom_line(color="blue", size=3)+ 
  theme_classic() +
  labs(x = "Temperature (C)", y = "Survivorship +/- SE") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size=16),
        strip.text = element_text(size=16))+
  geom_ribbon(data=Pachyheat, 
              alpha = 0.3, 
              aes(ymin = lwr, ymax = upr))+
  geom_point(data=Pachy.new,aes(x=TempC,y=Surv), size=2, alpha=0.5)

#Add to the master df
Heat.all.df <- rbind (Heat.all.df, Pachyheat)


#The 3 easier ones are done. Now to do Irrorata. Volker sent some code that should work

##Irrorata

Irrordat<-df.long%>%filter(Species=="Irrorata")
Irr.glm<-glm(Surv~TempC, data=Irrordat, family="binomial",control=glm.control(maxit=500))
summary(Irr.glm)
# issue is that we have complete /quasi complete separation, so there is only one temperature where we have 0 and 1 Temp2=97, everytihng <97 =1 and everything >97=0. So there is no variance to be explained
#Firth bias reduced log regression
library(logistf)
Irr.logistf<-logistf(formula=Surv~TempC, data=Irrordat, control=logistf.control(maxstep=10, maxit=500), pl=TRUE)
summary(Irr.logistf)

ndata <- with(Irrordat,data_frame(TempC = seq(min(TempC), max(TempC),length = 100)))
ndata<-add_column(ndata,fit = predict(Irr.logistf, newdata=ndata, type = "response"))

ggplot(ndata, aes(x=TempC, y=fit))+
  geom_vline(xintercept=35, linetype="dashed",size=2, color="grey")+
  geom_line(size=2, color="blue")+
  theme_classic()+
  geom_point(position=position_jitter(width=0.05,height=0.0), data=Irrordat, aes(x=TempC, y=Surv))

#That worked. Gonna need to add some columns to this data set in order to add it to the master df

Heat.all.df <- ndata%>%
  mutate(fit_link = NA,
         se_link = NA,
         fit_resp = NA,
         upr = fit,
         lwr = fit,
         Species = "Irrorata")%>%
  rbind(Heat.all.df)




#Calculate LT50 and add to a separate df
library(MASS)

#Trying to get LD50 standard errors by manually bootstrapping

numsim <- 10 #number of simulations to run

#dataframe to insert the subsampled observations

Buen.rand <- data.frame(Trial = NA,
                        Chamber = NA,
                        Bath = NA,
                        TempC = NA,
                        Species = NA,
                        Totalnum = NA,
                        Surv = NA,
                        sim = NA) 

#Function to run the bootstrap 

getLT <- function(species, numsim){    
  
  df <- data.frame(LT50 = NULL)  
  
  for(i in 1:numsim){
    
    n <- df.long%>%
      filter(Species == species)%>%
      nrow()
    
    rand.df <- df.long%>%
      filter(Species == species)%>% 
      sample_n(size = n, replace = T)%>%
      mutate(sim = i)
    
    model <- rand.df%>%
      filter(sim == i)%>%
      glm(data = ., Surv ~ TempC, family = "binomial")
    
    df <- data.frame(LT50 = unname(dose.p(model))[1],
                     Species = species)%>%
      rbind(df)
    
    
    
  }
  return(df)
}

#Run the function for each

LTdf <- data.frame(rbind(getLT("Buenoa", 1000),
                         getLT("Indica", 1000),
                         #getLT("Irrorata", 1000), #Irrorata uses a different function to model, so I'll do that separately below
                         getLT("Copto", 1000),
                         getLT("Pachy", 1000)))

#Bootstrap for Irrorata

Irr.LT <- for(i in 1:10){
  
  n <- df.long%>%
    filter(Species == "Irrorata")%>%
    nrow()
  
  rand.df <- df.long%>%
    filter(Species == "Irrorata")%>% 
    sample_n(size = n, replace = T)%>%
    mutate(sim = i)
  
  model <- rand.df%>%
    filter(sim == i)%>%
    logistf(formula = Surv ~ TempC, control = logistf.control(maxstep=10, maxit=500), pl=TRUE)
  
  df <- data.frame(LT50 = unname(dose.p(model))[1],
                   Species = "Irrorata")%>%
    rbind(df)
  
  df
}

LT.avg.df <- LTdf%>%
  group_by(Species)%>%
  dplyr::summarise(TempC = mean(LT50),
                   lwrCI = quantile(LT50, 0.025),
                   uprCI = quantile(LT50, 0.975))



#Graph all species

##Change Species in Heat.all.df to factor and relevel for graphing
Heat.all.df <- mutate_at(Heat.all.df, vars(Species), factor)%>%
  mutate(Species = fct_relevel(Species, c("Buenoa", "Indica", "Irrorata", "Copto", "Pachy")))

#Colorblind friendly palette
cbPalette <- c("#CC79A7", "#78C1EA", "#009E73", "#E69F00", "#D55E00", "#0072B2")

Heat.all.df%>%
  drop_na(Species)%>%
  ggplot(aes(x = TempC, y = fit, color = Species, fill = Species)) +
  geom_line(size=1) +
  facet_wrap(~Species, labeller = labeller(Species = c("Buenoa" = "Buenoa sp.", 
                                                       "Indica" = "N. indica", 
                                                       "Irrorata" = "N. irrorata", 
                                                       "Copto" = "C. loticus", 
                                                       "Pachy" = "P. longipennis",
                                                       "Tramea" = "T. carolina"))) +
  geom_point(data = df.long, aes(y = Surv), size = 2) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3, linetype = 0) +
  theme_classic() +
  labs(x = "Temperature (C)", y = "Survivorship +/- SE") +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size=18),
        strip.text = element_text(size=16, face = "italic"),
        legend.position = 0) +
  #geom_vline(xintercept = LTdf$TempC, linetype = "dotted", size = 2, alpha = 0.7) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  geom_point(data = LT.avg.df, aes(x = TempC, y = .5), color = "black", size = 2) + 
  geom_errorbarh(data = LT.avg.df, aes(xmin = lwrCI, xmax = uprCI, y = .5), color = "black", height = 0.1) 


ggsave("Figures/Heat.tolerance.pdf", width = 13.32, height = 7.27)




#####################################


#Here's a couple other methods for calculating the LT50. I might want to go with one of these

##Trying a function that calculates LD50 using a logistic regression


library(HelpersMG)

#First I'll try with Buenoa
df%>%
  filter(Species == "Buenoa")%>%
  data.frame()%>%
  mutate(alive = Numsurv, dead = Numdead, N = Totalnum, doses = TempC)%>%
  LD50()



#Trying a different function too, just to compare
library(ecotox)


df.long%>%
  filter(Species == "Buenoa")%>%
  LC_logit(Surv ~ TempC, data = ., p = 50)


#Trying the dose.p function in Mass
##Write a function to create a data frame of all species

getLT <- function(species, model){
  data.frame(Species = species,
             TempC = unname(dose.p(model))[1],
             SE = dose.p(model)%>%attr("SE")%>%unname())
  
}

###Run function for each species and crate df.

LTdf <- rbind(getLT("Buenoa", Buen.glm),
              getLT("Indica", Indica.glm),
              getLT("Irrorata", Irr.logistf),
              getLT("Copto", Copto.glm),
              getLT("Pachy", Pachy.glm))%>%
  mutate(Species = fct_relevel(Species, c("Buenoa", "Indica", "Irrorata", "Copto", "Pachy"))) #Relevel for graphing


