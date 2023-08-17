
library(tidyverse)

#load the raw data, add a variable of the number of dead daphnia
feedRate.df <- read.csv("Data/pred.hier.feed.csv", header=T)%>%
  mutate(dead = 100 - surv)%>%
  mutate_at(vars(c(trial, bath, location)), funs(factor))


#Replace surv > 100 with 100. Any values above 100 were from counting error
feedRate.df$surv[feedRate.df$surv > 100]<-100


#Next I'll do the actual GLM to get the predicted survivor values for the controls.  Actually, not that I think about it a GLM wouldn't work because I expect the relationship to be nonlinear.  A GAM might be better, though this may be dependent on how I ultimately analyze the feeding rates.  I'll stick with GAM for now and revisit this down the line.

#Add column for number dead for the binomial glm's

feedRate.df <- feedRate.df%>%
  mutate(dead = 100 - surv)


#Calculate background mortality with binomial glmer

library(lme4)

#Filter all data except control treatments
cont.df <- feedRate.df%>%
  filter(pred == "Control")


cont.glmer <- cont.df%>%
  glmer(data = ., cbind(surv, dead) ~ temp  + (1|trial) + (1|bath) + (1|location), family = "binomial")

plot(cont.glmer)


#Plot with data
cont.df%>%
  mutate(predicted = predict(cont.glmer, newdata = cont.df, 
                             type = "response", re.form = ~0))%>%
  ggplot(aes(x = temp, y = surv)) +
  geom_line(aes(y = predicted * 100)) +
  geom_point()

#Add predicted values to control df
cont.df <- cont.df%>%
  mutate(backmortProp = 1 - predict(cont.glmer, newdata = cont.df,
                                    type = "response", re.form = ~0),
         backmortNum = backmortProp * 100)


#Subtract these values from 100 to get predicted background mortality for each of these temperatures
feed.data <-cont.df%>%
  dplyr::select(trial, bath, backmortProp, backmortNum)%>%
  right_join(feedRate.df, by = c("trial","bath"), keep = F)%>%
  mutate(Numeaten = dead-backmortNum)

##Convert negative feeding rate values to 0's.  These are the result of counting errors
feed.data$Numeaten <- (feed.data$Numeaten + abs(feed.data$Numeaten))/2

##Remove cases where the predators didn't eat anything
feed.data <- filter(feed.data, Numeaten != 0)


#Analyze the feeding rates

##Create a function to select the best model using AICc

library(MuMIn)
library(lmtest)
selectmodel<- function(species){
  df <- feed.data%>%
    drop_na(Numeaten)%>%
    filter(pred == species)
  
  
  
  model1 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2)  + scale(predmass) + (1|trial) + (1|bath) + (1|location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model2 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2)  + scale(predmass) + (1|trial) + (1|bath), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model3 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2)  + scale(predmass) + (1|trial), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model4 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2)  + scale(predmass) + (1|trial)  + (1|location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model5 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2)  + scale(predmass) + (1|location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model6 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2)  + scale(predmass) + (1|bath), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model7 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + (1|trial) + (1|bath) + (1|location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model8 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + (1|trial) + (1|bath), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model9 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + (1|trial), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model10 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + (1|trial)  + (1|location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model11 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + (1|location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model12 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + (1|bath), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model13 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp  + scale(predmass) + (1|trial) + (1|bath) + (1|location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model14 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp  + scale(predmass) + (1|trial) + (1|bath), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model15 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp  + scale(predmass) + (1|trial), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model16 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp  + scale(predmass) + (1|trial)  + (1|location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model17 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp  + scale(predmass) + (1|location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model18 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp  + scale(predmass) + (1|bath), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model19 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp + (1|trial) + (1|bath) + (1|location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model20 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp + (1|trial) + (1|bath), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model21 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp + (1|trial), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model22 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp + (1|trial)  + (1|location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model23 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp + (1|location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  model24 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp + (1|bath), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))
  
  
  
  
  AICdf <- data.frame(model = c("model1", "model2", "model3", "model4", "model5", "model6", 
                                "model7", "model8", "model9", "model10", "model11", "model12",
                                "model13", "model14", "model15", "model16", "model17", "model18", 
                                "model19", "model20", "model21", "model22", "model23", "model24"),
                      AIC = c(AICc(model1), AICc(model2), AICc(model3), AICc(model4), AICc(model5), AICc(model6),
                              AICc(model7), AICc(model8), AICc(model9), AICc(model10), AICc(model11), AICc(model12), 
                              AICc(model13), AICc(model14), AICc(model15), AICc(model16), AICc(model17), AICc(model18),
                              AICc(model19), AICc(model20), AICc(model21), AICc(model22), AICc(model23), AICc(model24)))%>%
    arrange(AIC)%>%
    mutate(delta = AIC - .[1,2])%>%
    return()
  
  list(AICdf,
       lrtest(get(AICdf[1,1]), get(AICdf[2,1])),
       lrtest(get(AICdf[1,1]), get(AICdf[3,1])),
       lrtest(get(AICdf[1,1]), get(AICdf[4,1])))
  
}


BuenAIC <- selectmodel("Buenoa")
IndAIC <- selectmodel("Indica")
IrrAIC <- selectmodel("Irrorata")
CoptoAIC <- selectmodel("Copto")
PachyAIC <- selectmodel("Pachy")
TramAIC <- selectmodel("Tramea")

#Run the best models for each species

Buen.glmer <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(pred == "Buenoa")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp  + scale(predmass) + (1|trial) + (1|bath), family = "binomial", data = ., control = glmerControl(optimizer = "bobyqa"))

Ind.glmer <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(pred == "Indica")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + scale(predmass) + (1|bath) + (1|trial) + (1|location), family = "binomial", 
        data = ., control = glmerControl(optimizer = "bobyqa"))

Irr.glmer <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(pred == "Irrorata")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + (1|trial) + (1|bath) + (1|location), family = "binomial", 
        data = ., control = glmerControl(optimizer = "bobyqa"))

Copto.glmer <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(pred == "Copto")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2)  + (1|trial) + (1|bath) + (1|location), family = "binomial", 
        data = ., control = glmerControl(optimizer = "bobyqa"))

Pachy.glmer <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(pred == "Pachy")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + (1|trial) + (1|bath) + (1|location), family = "binomial", 
        data = ., control = glmerControl(optimizer = "bobyqa"))

Tram.glmer <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(pred == "Tramea")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + (1|bath), family = "binomial", 
        data = ., control = glmerControl(optimizer = "bobyqa"))

#Add predicted values and plot them

library(AICcmodavg) #Package for generating standard errors from glmer objects

#library(merTools) #Package with function for calculating prediction intervals from glmer.
##Create function for generating predicted values

getpredict <- function(model, species, subdata){
  subdata%>%
    #cbind(predictInterval(model, newdata = subdata, type = "probability", 
                          #which = "fixed", include.resid.var = F))%>% #Uncomment if prediction interval desired instead of se
    cbind(predictSE.mer(model, se.fit = T, newdata = subdata, type = "response"))%>% #Generate standard errors
    mutate(lwr = fit - se.fit, upr = fit + se.fit)%>%
    dplyr::select(-se.fit)%>%
    return()
  
}


##Buenoa
###Create the subset data
Buendata <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(pred == "Buenoa")

###Set predmass = to mean value
Buendata$predmass <- mean(Buendata$predmass, na.rm = T)

###Calculate predicted values
Buen.predict <- getpredict(Buen.glmer, "Buenoa", Buendata)


##Indica
###Create the subset of data
Indata <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(pred == "Indica")

###Set predmass = to mean value
Indata$predmass <- mean(Indata$predmass, na.rm = T)


###Calculate the predicted values
Ind.predict <- getpredict(Ind.glmer, "Indica", Indata)

##Irrorata
###Create the subset of data
Irrdata <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(pred == "Irrorata")

###Set predmass = to mean value
Irrdata$predmass <- mean(Irrdata$predmass, na.rm = T)

###Calculate the predicted values
Irr.predict <- getpredict(Irr.glmer, "Irrorata", Irrdata)


#Copto
###Create the subset of data
Coptodata <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(pred == "Copto")

###Set predmass = to mean value
Coptodata$predmass <- mean(Coptodata$predmass, na.rm = T)

###Calculate the predicted values
Copto.predict <- getpredict(Copto.glmer, "Copto", Coptodata)

#Pachy
###Create the subset of data
Pachydata <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(pred == "Pachy")

###Set predmass = to mean value
Pachydata$predmass <- mean(Pachydata$predmass, na.rm = T)

###Calculate the predicted values
Pachy.predict <- getpredict(Pachy.glmer, "Pachy", Pachydata)

#Tramea
###Create the subset of data
Tramdata <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(pred == "Tramea")

###Set predmass = to mean value
Tramdata$predmass <- mean(Tramdata$predmass, na.rm = T)

###Calculate the predicted values
Tram.predict <- getpredict(Tram.glmer, "Tramea", Tramdata)

##Combine each of these to a final data set
final.data <- rbind(Buen.predict, Irr.predict, Ind.predict, Copto.predict, Pachy.predict, Tram.predict)%>%
  mutate(pred = factor(pred))%>%
  mutate(pred = fct_relevel(pred, c("Buenoa", "Indica", "Irrorata", "Copto", "Pachy", "Tramea")))

###################

#Now graph it


##colorblind friendly palette

cbPalette <- c("#CC79A7", "#78C1EA", "#009E73", "#E69F00", "#D55E00", "#0072B2")

##Species labels
spLabs <- c("Buenoa" = "B. scimitra", 
            "Indica" = "N. indica", 
            "Irrorata" = "N. irrorata", 
            "Copto" = "C. loticus", 
            "Pachy" = "P. longipennis",
            "Tramea" = "T. carolina")

###Faceted by species
final.data%>%
  ggplot(aes(x = temp, y = fit*100, color = pred, fill = pred)) +
  facet_wrap(~pred, labeller = labeller(pred = spLabs)) +
  geom_point(aes(y = Numeaten)) +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = 0.25, aes(ymin = lwr*100, ymax = upr*100), linetype = 0) +
  theme_classic() +
  labs(x = "Temperature (C)", y = "Number of prey eaten ± SE") +   
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size=24),
        strip.text = element_text(size=22, face = "italic"),
        legend.position = 0) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette)

#Uncomment to save
ggsave("Figures/Feedresults.jpeg", width = 13.32, height = 7.27)

###All species in one panel WITHOUT error bars

#Using color so here's a colorblind friendly palette
final.data%>%
  ggplot(aes(x = temp, y = fit*100)) +
  #geom_point(aes(y = Numeaten)) +
  geom_line(aes(color = pred), size = 2) +
  theme_classic() +
  labs(x = "Temperature (C)", y = "Number of prey eaten ± SE") +   
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size=18),
        legend.text = element_text(size = 14, face = "italic"),
        legend.title = element_text(size = 16),
        legend.position = "none") +
  scale_color_manual(values = cbPalette, name = "species",
                     labels = spLabs) +
  scale_fill_manual(values = cbPalette) 

#Uncomment to save
#ggsave("Figures/Feedresults2.jpeg", width = 6.5, height = 6.0)


#All species in one panel WITH error bars

final.data%>%
  ggplot(aes(x = temp, y = fit*100)) +
  geom_line(aes(color = pred), size = 2) +
  geom_ribbon(alpha = 0.25, aes(ymin = lwr*100, ymax = upr*100, fill = pred), show.legend = F) +
  theme_classic() +
  labs(x = "Temperature (C)", y = "Number of prey eaten ± SE") +   
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size=18),
        legend.text = element_text(size = 14, face = "italic"),
        legend.title = element_text(size = 16)) +
  scale_color_manual(values = cbPalette, name = "species",
                     labels = spLabs) +
  scale_fill_manual(values = cbPalette) 

#Uncomment to save
#ggsave("Figures/Feedresults3.jpeg", width = 9.5, height = 5.9)
















