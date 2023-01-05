

#First, a small bit of data cleaning.  I need to make temp2 values equal to temp1 values for all rows from trial 5 plus baths 7 and 8 and trial 8

library(tidyverse)
#load the raw data
surv.df <- read.csv("Data/Data_for_r1.csv",header=T)%>%
  mutate(Dead = 100 - Surv,
         Temp = (Temp2 - 32)*5/9)%>% #This is the temperature I'll use for analyses
  mutate_at(vars(c(Trial, Bath, Location)), funs(factor))

#Load predator body mass data
mass.df <- read.csv("Pred\ Dry\ Mass.csv", header = T)

#Check predator mass df for duplicates
mass.df%>%
  group_by(Trial, Bath, Location)%>%
  dplyr::mutate(dupe = n()>1)%>%
  dplyr::filter(dupe == T)

mass.df%>%
  group_by(Trial, Bath, Location)%>%
  filter(n()<2)

#add the predator mass data
surv.mass.df <- mass.df%>%
  group_by(Trial, Bath, Location)%>%
  filter(n() <2)%>% #Remove instances where there are two rows with different weights for the same predator
  dplyr::select(c(Trial, Bath, Location, Predmass))%>%
  mutate_at(vars(c(Trial, Bath, Location)), funs(factor))%>%
  right_join(surv.df, by = c("Trial", "Bath", "Location"))
  

#Replace NA values in Predmass with species averages

surv.mass.df <- surv.mass.df%>%
  group_by(Pred)%>%
  mutate(Predmass = replace_na(Predmass, mean(Predmass, na.rm = T)))

#I need to replace temp2 with temp1 for trial 5 because I didn't record it. For now I replaced it in the excel file, but I might want to do it in the code 

#Next I'll do the actual GLM to get the predicted survivor values for the controls.  Actually, not that I think about it a GLM wouldn't work because I expect the relationship to be nonlinear.  A GAM might be better, though this may be dependent on how I ultimately analyze the feeding rates.  I'll stick with GAM for now and revisit this down the line.

#Replace Surv > 100 with 100
surv.mass.df$Surv[surv.mass.df$Surv > 100]<-100

#Add column for number dead for the binomial glm's

surv.mass.df <- surv.mass.df%>%
  mutate(Dead = 100 - Surv)


#Calculate background mortality with binomial glmer

library(lme4)

#Filter all data except control treatments
cont.df <- surv.mass.df%>%
  filter(Pred == "Control")


cont.glmer <- cont.df%>%
  glmer(data = ., cbind(Surv, Dead) ~ Temp  + (1|Trial) + (1|Bath) + (1|Location), family = "binomial")

plot(cont.glmer)


#Plot with data
cont.df%>%
  mutate(predicted = predict(cont.glmer, newdata = cont.df, 
                             type = "response", re.form = ~0))%>%
  ggplot(aes(x = Temp, y = Surv)) +
  geom_line(aes(y = predicted * 100)) +
  geom_point()

#Add predicted values to control df
cont.df <- cont.df%>%
  mutate(backmortProp = 1 - predict(cont.glmer, newdata = cont.df,
                                    type = "response", re.form = ~0),
         backmortNum = backmortProp * 100)


#Next, I'll need to subtract these values from 100 to get predicted background mortality for each of these temperatures.  Ultimately these background mortality values will be subtracted from the mortality values observed in the predator treatments.  I think the best way to do this is to attach the background mortality values to the data set.

feed.data <-cont.df%>%
  dplyr::select(Trial, Bath, backmortProp, backmortNum)%>%
  right_join(surv.mass.df, by = c("Trial","Bath"), keep = F)%>%
  mutate(Numeaten = Dead-backmortNum)
  
  

##Convert negative feeding rate values to 0's.  These are the result of counting errors
feed.data$Numeaten <- (feed.data$Numeaten + abs(feed.data$Numeaten))/2

##Remove cases where the predators didn't eat anything
feed.data <- filter(feed.data, Numeaten != 0)


#I have the final data set of number eaten minus the predicted background mortality counts. Now to run glmer's on them. I'll do a separate model per species. I'll compare single and double term models to look for nonlinearity


##Trying out a global model with all predators and species as fixed effect

All.glmer <- feed.data%>%
  filter(Pred != "Control")%>%
  drop_na(Numeaten)%>%
  glmer(cbind(round(Numeaten), round(100 - Numeaten)) ~ Temp + Pred + I(Temp^2) + Pred*Temp + Pred*I(Temp^2) +
                (1|Predmass), family = "binomial", data = .,
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

#Failed to converge

#Trying glm
All.glm <- feed.data%>%
  filter(Pred != "Control")%>%
  na.omit()%>%
  glm(cbind(round(Numeaten), round(100 - Numeaten)) ~ Pred:Temp + Pred:I(Temp^2) +  Pred + Pred:Predmass, family = "binomial", data = ., na.action = "na.pass")
      
library(MuMIn)
#Check AIC values for different models
Model1 = feed.data%>%
  filter(Pred != "Control")%>%
  drop_na(Numeaten)%>%
  glm(cbind(round(Numeaten), round(100 - Numeaten)) ~ Temp + I(Temp^2) +  Pred + Predmass, family = "binomial", data = .)
Model2 = feed.data%>%
  filter(Pred != "Control")%>%
  drop_na(Numeaten)%>%
  glm(cbind(round(Numeaten), round(100 - Numeaten)) ~ Pred:Temp + I(Temp^2) +  Pred + Predmass, family = "binomial", data = .)
Model3 = feed.data%>%
  filter(Pred != "Control")%>%
  drop_na(Numeaten)%>%
  glm(cbind(round(Numeaten), round(100 - Numeaten)) ~ Pred:Temp + Pred:I(Temp^2) +  Pred + Predmass, family = "binomial", data = .)

#This one has lowest AIC
Model4 = feed.data%>%
  filter(Pred != "Control")%>%
  drop_na(Numeaten)%>%
  glm(cbind(round(Numeaten), round(100 - Numeaten)) ~ Pred:Temp + Pred:I(Temp^2) +  Pred + Pred:Predmass, family = "binomial", data = .)


Model5 = feed.data%>%
  filter(Pred != "Control")%>%
  drop_na(Numeaten)%>%
  glm(cbind(round(Numeaten), round(100 - Numeaten)) ~ Pred:Temp + I(Temp^2) +  Pred + Pred:Predmass, family = "binomial", data = .)
Model6 = feed.data%>%
  filter(Pred != "Control")%>%
  drop_na(Numeaten)%>%
  glm(cbind(round(Numeaten), round(100 - Numeaten)) ~ Temp + I(Temp^2) +  Pred + Pred:Predmass, family = "binomial", data = .)
Model7 = feed.data%>%
  filter(Pred != "Control")%>%
  drop_na(Numeaten)%>%
  glm(cbind(round(Numeaten), round(100 - Numeaten)) ~ Predmass:Temp + Predmass:I(Temp^2) +  Pred + Predmass, family = "binomial", data = .)
Model8 = feed.data%>%
  filter(Pred != "Control")%>%
  drop_na(Numeaten)%>%
  glm(cbind(round(Numeaten), round(100 - Numeaten)) ~ Predmass:Temp + Predmass:I(Temp^2) +  Pred + Predmass, family = "binomial", data = .)

Model9 = feed.data%>%
  filter(Pred != "Control")%>%
  drop_na(Numeaten)%>%
  glm(cbind(round(Numeaten), round(100 - Numeaten)) ~ Predmass:Temp + Predmass:I(Temp^2) +  Pred + Pred:Predmass + Predmass, family = "binomial", data = .)
Model10 = feed.data%>%
  filter(Pred != "Control")%>%
  drop_na(Numeaten)%>%
  glm(cbind(round(Numeaten), round(100 - Numeaten)) ~ Pred:Predmass:Temp + I(Temp^2) +  Pred + Predmass + Predmass:Pred, family = "binomial", data = .)
Model11 = feed.data%>%
  filter(Pred != "Control")%>%
  drop_na(Numeaten)%>%
  glm(cbind(round(Numeaten), round(100 - Numeaten)) ~ Pred:Predmass:Temp + I(Temp^2) +  Pred + Predmass, family = "binomial", data = .)
Model12 = feed.data%>%
  filter(Pred != "Control")%>%
  drop_na(Numeaten)%>%
  glm(cbind(round(Numeaten), round(100 - Numeaten)) ~ Pred:Predmass:Temp + Pred:Predmass:I(Temp^2) +  Pred + Predmass, family = "binomial", data = .)
Model13 = feed.data%>%
  filter(Pred != "Control")%>%
  drop_na(Numeaten)%>%
  glm(cbind(round(Numeaten), round(100 - Numeaten)) ~ Pred:Predmass:Temp + Pred:Predmass:I(Temp^2) +  Pred + Predmass + Pred:Predmass, family = "binomial", data = .)

Model14 = feed.data%>%
  filter(Pred != "Control")%>%
  drop_na(Numeaten)%>%
  glm(cbind(round(Numeaten), round(100 - Numeaten)) ~ Pred:Predmass:Temp + Pred:Predmass:I(Temp^2) +  Pred + Predmass:Pred, family = "binomial", data = .)


feed.data%>%
  drop_na(Numeaten)%>%
  filter(Pred == "Irrorata")%>%
  ggplot(aes(x = Temp, y = Numeaten)) +
           #geom_point(position = "jitter") +
           #geom_jitter(width = 1) +
           geom_point(aes(y = predict(All.glm, newdata = drop_na(filter(feed.data, Pred == "Irrorata"), Numeaten), type = "response")*100))
 
feed.data%>%
  drop_na(Numeaten)%>%
  filter(Pred == "Irrorata")
  plot(x = filter(drop_na(feed.data), Pred == "Irrorata")$Temp, y = predict(All.glm, newdata = drop_na(filter(feed.data, Pred == "Irrorata"))))

  
#############

#The global model isn't converging. Here's the species-specific models

##Create a function for the glmer's. The function takes the predator and number of terms for polynomial to use


runglmer <- function(predator, terms){
  feed.data%>%
    drop_na(Numeaten)%>%
    filter(Pred == predator)%>%
    glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, terms)  + scale(Predmass) + (1|Trial) + (1|Bath), family = "binomial", data = ., control = glmerControl(optimizer = "bobyqa"))%>%
    return()
  
  
}


#Going to try using the best model for each species rather than keeping them all the same. Here I'm just testing some code to find the best way to do this across species

library(MuMIn)

selectmodel<- function(Species){
  df <- feed.data%>%
    drop_na(Numeaten)%>%
    filter(Pred == Species)
    
  

model1 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, 2)  + scale(Predmass) + (1|Trial) + (1|Bath) + (1|Location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model2 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, 2)  + scale(Predmass) + (1|Trial) + (1|Bath), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model3 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, 2)  + scale(Predmass) + (1|Trial), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model4 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, 2)  + scale(Predmass) + (1|Trial)  + (1|Location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model5 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, 2)  + scale(Predmass) + (1|Location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model6 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, 2)  + scale(Predmass) + (1|Bath), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model7 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, 2) + (1|Trial) + (1|Bath) + (1|Location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model8 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, 2) + (1|Trial) + (1|Bath), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model9 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, 2) + (1|Trial), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model10 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, 2) + (1|Trial)  + (1|Location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model11 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, 2) + (1|Location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model12 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, 2) + (1|Bath), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model13 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp  + scale(Predmass) + (1|Trial) + (1|Bath) + (1|Location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model14 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp  + scale(Predmass) + (1|Trial) + (1|Bath), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model15 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp  + scale(Predmass) + (1|Trial), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model16 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp  + scale(Predmass) + (1|Trial)  + (1|Location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model17 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp  + scale(Predmass) + (1|Location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model18 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp  + scale(Predmass) + (1|Bath), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model19 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp + (1|Trial) + (1|Bath) + (1|Location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model20 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp + (1|Trial) + (1|Bath), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model21 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp + (1|Trial), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model22 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp + (1|Trial)  + (1|Location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model23 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp + (1|Location), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))

model24 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp + (1|Bath), family = "binomial", data = df, control = glmerControl(optimizer = "bobyqa"))




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
  filter(Pred == "Buenoa")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp  + scale(Predmass) + (1|Trial) + (1|Bath), family = "binomial", data = ., control = glmerControl(optimizer = "bobyqa"))

Ind.glmer <- feed.data%>%
   drop_na(Numeaten)%>%
  filter(Pred == "Indica")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, 2) + scale(Predmass) + (1|Bath) + (1|Trial) + (1|Location), family = "binomial", 
        data = ., control = glmerControl(optimizer = "bobyqa"))

Irr.glmer <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(Pred == "Irrorata")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, 2) + (1|Trial) + (1|Bath) + (1|Location), family = "binomial", 
        data = ., control = glmerControl(optimizer = "bobyqa"))

Copto.glmer <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(Pred == "Copto")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, 2)  + (1|Trial) + (1|Bath) + (1|Location), family = "binomial", 
        data = ., control = glmerControl(optimizer = "bobyqa"))
  
Pachy.glmer <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(Pred == "Pachy")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, 2) + (1|Trial) + (1|Bath) + (1|Location), family = "binomial", 
        data = ., control = glmerControl(optimizer = "bobyqa"))

Tram.glmer <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(Pred == "Tramea")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, 2) + (1|Bath), family = "binomial", 
        data = ., control = glmerControl(optimizer = "bobyqa"))



#Add predicted values and plot them. I'll subset the data set into df for each species, calculate the predicted values, then join those to the final data frame for graphing

library(merTools) #Package with function for calculating predicted SE from glmer

feed.data <- as.data.frame(feed.data)


##Create function for generating predicted values

getpredict <- function(model, species, subdata){
  feed.data%>%
    drop_na(Numeaten)%>%
    filter(Pred == species)%>%
    cbind(predictInterval(model, newdata = subdata, type = "probability", 
                          which = "fixed", include.resid.var = F))%>%
    return()
  
}


##Buenoa
###Create the subset data
Buendata <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(Pred == "Buenoa")

###Set predmass = to mean value
Buendata$Predmass <- mean(Buendata$Predmass, na.rm = T)

###Calculate predicted values
Buen.predict <- getpredict(Buen.glmer, "Buenoa", Buendata)


##Indica
###Create the subset of data
Indata <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(Pred == "Indica")

###Set predmass = to mean value
Indata$Predmass <- mean(Indata$Predmass, na.rm = T)


###Calculate the predicted values
Ind.predict <- getpredict(Ind.glmer, "Indica", Indata)

##Irrorata
###Create the subset of data
Irrdata <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(Pred == "Irrorata")

###Set predmass = to mean value
Irrdata$Predmass <- mean(Irrdata$Predmass, na.rm = T)

###Calculate the predicted values
Irr.predict <- getpredict(Irr.glmer, "Irrorata", Irrdata)


#Copto
###Create the subset of data
Coptodata <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(Pred == "Copto")

###Set predmass = to mean value
Coptodata$Predmass <- mean(Coptodata$Predmass, na.rm = T)

###Calculate the predicted values
Copto.predict <- getpredict(Copto.glmer, "Copto", Coptodata)

#Pachy
###Create the subset of data
Pachydata <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(Pred == "Pachy")

###Set predmass = to mean value
Pachydata$Predmass <- mean(Pachydata$Predmass, na.rm = T)

###Calculate the predicted values
Pachy.predict <- getpredict(Pachy.glmer, "Pachy", Pachydata)

#Tramea
###Create the subset of data
Tramdata <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(Pred == "Tramea")

###Set predmass = to mean value
Tramdata$Predmass <- mean(Tramdata$Predmass, na.rm = T)

###Calculate the predicted values
Tram.predict <- getpredict(Tram.glmer, "Tramea", Tramdata)

##Combine each of these to a final data set
final.data <- rbind(Buen.predict, Irr.predict, Ind.predict, Copto.predict, Pachy.predict, Tram.predict)%>%
  mutate(Pred = factor(Pred))%>%
  mutate(Pred = fct_relevel(Pred, c("Buenoa", "Indica", "Irrorata", "Copto", "Pachy", "Tramea")))

###################

#Now graph it
  

##colorblind friendly palette

cbPalette <- c("#CC79A7", "#78C1EA", "#009E73", "#E69F00", "#D55E00", "#0072B2")

###Faceted by species
final.data%>%
  ggplot(aes(x = Temp, y = fit*100, color = Pred, fill = Pred)) +
  facet_wrap(~Pred, labeller = labeller(Pred = c("Buenoa" = "Buenoa sp.", 
                                                 "Indica" = "N. indica", 
                                                 "Irrorata" = "N. irrorata", 
                                                 "Copto" = "C. loticus", 
                                                 "Pachy" = "P. longipennis",
                                                 "Tramea" = "T. carolina"))) +
  geom_point(aes(y = Numeaten)) +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = 0.25, aes(ymin = lwr*100, ymax = upr*100), linetype = 0) +
  theme_classic() +
  labs(x = "Temperature (C)", y = "Number of Prey Eaten +/- CI") +   
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size=18),
        strip.text = element_text(size=16, face = "italic"),
        legend.position = 0) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette)


ggsave("Figures/Feedresults.pdf", width = 13.32, height = 7.27)

###All species in one panel WITHOUT error bars

#Using color so here's a colorblind friendly palette


final.data%>%
  ggplot(aes(x = Temp, y = fit*100)) +
  #geom_point(aes(y = Numeaten)) +
  geom_line(aes(color = Pred), size = 2) +
  theme_classic() +
  labs(x = "Temperature (C)", y = "# Eaten") +   
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size=18),
        legend.text = element_text(size = 14, face = "italic"),
        legend.title = element_text(size = 16),
        legend.position = "none") +
  scale_color_manual(values = cbPalette, name = "Species",
                     labels = c("Buenoa" = "Buenoa sp.", 
                                "Indica" = "N. indica", 
                                "Irrorata" = "N. irrorata", 
                                "Copto" = "C. loticus", 
                                "Pachy" = "P. longipennis",
                                "Tramea" = "T. carolina")) +
  scale_fill_manual(values = cbPalette) 


ggsave("Figures/Feedresults2.pdf", width = 6.5, height = 6.0)


#All species in one panel WITH error bars

final.data%>%
  ggplot(aes(x = Temp, y = fit*100)) +
  geom_line(aes(color = Pred), size = 2) +
  geom_ribbon(alpha = 0.25, aes(ymin = lwr*100, ymax = upr*100, fill = Pred), show.legend = F) +
  theme_classic() +
  labs(x = "Temperature (C)", y = "# Eaten +/- SE") +   
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size=18),
        legend.text = element_text(size = 14, face = "italic"),
        legend.title = element_text(size = 16)) +
  scale_color_manual(values = cbPalette, name = "Species",
                     labels = c("Buenoa" = "Buenoa sp.", 
                                "Indica" = "N. indica", 
                                "Irrorata" = "N. irrorata", 
                                "Copto" = "C. loticus", 
                                "Pachy" = "P. longipennis",
                                "Tramea" = "T. carolina")) +
  scale_fill_manual(values = cbPalette) 

ggsave("Figures/Feedresults3.pdf", width = 9.5, height = 5.9)
















#############################################
#Old code, just keeping here just in case.

##I'll try GLM's too to compare the models
runglm <- function(species, terms){
  feed.data%>%
    na.omit()%>%
    filter(Pred == species)%>%
    glm(data = .,cbind(round(Numeaten), round(100 - Numeaten)) ~ poly(Temp, terms) + Predmass,
        family = "binomial")%>%
    return()
}




##Looking for effects of each random effect
feed.data%>%
  drop_na(Numeaten)%>%
  filter(Pred == "Buenoa")%>%
  ggplot(aes(x = Temp, y = Numeaten, color = Location)) + 
  #geom_smooth() +
  geom_point() +
  geom_jitter(width = 1)

library(car)
##Testing for significance in random effects by including them as fixed. Just gonna run these a bunch of times with different predators and different effects as fixed
feed.data%>%
  drop_na(Numeaten)%>%
  filter(Pred == "Pachy")%>%
  glmer(data = ., cbind(round(Numeaten), round(100 - Numeaten)) ~ poly(Temp, 2) +
          (1|Trial) + (1|Bath) + (1|Location), family = "binomial")%>%
  Anova(type = 3)

##They're all significant for all species. I'm moving forward with the full models


#Run the 1 and 2 term models for each predator




##Buenoa - first order glmer. Delta AIC = 1.0154

Buen.glmer <- runglmer("Buenoa", 2)

data.frame(glmer1 = AIC(runglmer("Buenoa", 1)),
           glmer2 = AIC(runglmer("Buenoa", 2)),
           glm1 = AIC(runglm("Buenoa", 1)),
           glm2 = AIC(runglm("Buenoa", 2)))

##Indica - second order glmer. Delta AIC = 20.124

data.frame(glmer1 = AIC(runglmer("Indica", 1)),
           glmer2 = AIC(runglmer("Indica", 2)),
           glm1 = AIC(runglm("Indica", 1)),
           glm2 = AIC(runglm("Indica", 2)))

##Likelihood ratio test
library(lmtest)

lrtest(runglmer("Indica", 1), runglmer("Indica", 2)) #chi-sq p value sig

Ind.glmer <- runglmer("Indica", 2)

##Irrorata - second order glmer. Delta AIC = 38.8905
data.frame(glmer1 = AIC(runglmer("Irrorata", 1)),
           glmer2 = AIC(runglmer("Irrorata", 2)),
           glm1 = AIC(runglm("Irrorata", 1)),
           glm2 = AIC(runglm("Irrorata", 2)))

lrtest(runglmer("Irrorata", 1), runglmer("Irrorata", 2)) #p value sig

Irr.glmer <- runglmer("Irrorata", 2)

##Copto - second glmer. Delta AIC = 2.3712
data.frame(glmer1 = AIC(runglmer("Copto", 1)),
           glmer2 = AIC(runglmer("Copto", 2)),
           glm1 = AIC(runglm("Copto", 1)),
           glm2 = AIC(runglm("Copto", 2)))

lrtest(runglmer("Copto", 1), runglmer("Copto", 2)) #p value sig

Copto.glmer <- runglmer("Copto", 2)

##Pachy - first order glmer. Delta AIC = 0.2969
data.frame(glmer1 = AIC(runglmer("Pachy", 1)),
           glmer2 = AIC(runglmer("Pachy", 2)),
           glm1 = AIC(runglm("Pachy", 1)),
           glm2 = AIC(runglm("Pachy", 2)))

lrtest(runglmer("Pachy", 1), runglmer("Pachy", 2))# p value not sig

Pachy.glmer <- runglmer("Pachy", 2)


##Tramea 2nd order is over fit. Trying a simpler model for the polynomial

data.frame(glmer1 = AIC(runglmer("Tramea", 1)),
           glmer2 = AIC(runglmer("Tramea", 2)),
           glm1 = AIC(runglm("Tramea", 1)),
           glm2 = AIC(runglm("Tramea", 2)))
Tram.glmer <- runglmer("Tramea", 2)




#The Indica model is looking funky. It's underestimating in the middle. Here I'm just trying to work that out


#Add second order temperature values
Indata <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(Pred == "Indica")%>%
  mutate(Temp2 = Temp^2,
         roundNumeaten = round(Numeaten),
         roundNumsurv = round(100-Numeaten))



#Pasting the model here for tweaking
Ind.glmer.test <-feed.data%>%
  drop_na(Numeaten)%>%
  filter(Pred == "Indica")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp,2)  + scale(Predmass) + (1|Trial) + (1|Bath) + (1|Location), family = "binomial", data = ., control = glmerControl(optimizer = "bobyqa"))

#Get predicteds
Ind.predict.test <- getpredict(Ind.glmer.test, "Indica", Indata)

#Plot

feed.data%>%
  drop_na(Numeaten)%>%
  filter(Pred == "Indica")%>%
  cbind(predictInterval(Ind.glmer.test, newdata = ., type = "probability", 
                        which = "fixed", include.resid.var = F))%>%
  ggplot(aes(x = Temp, y = Numeaten)) +
  geom_point() +
  geom_line(aes(x = Temp, y = fit * 100))

averageObs(Ind.glmer.test)



Ind.predict%>%
  ggplot(aes(x = Bath, y = fit)) +
  geom_point()




#The data looks like it's a saturating curve. I'll try nls with a logistic growth formula

library(nlme)

# I'll need to include starter values for the parameters. I'll get estimates for r and a with coefficients from a linear model

lm(logit(Numeaten/50) ~ Temp, Indata)

#I'll use 0.03 for r and -5.7 for a. The data seems to saturate at around 50 eaten, so I'll use log(50) as the starting value for k 
#Run the nonlinear model

nls(Numeaten ~ I(k / 1 + exp(-r * (Temp - a))), data = Indata,
    start = list(k = log(50), r = 0.03, a = -5.7))

#I'm getting an error "singular gradient." Some people online suggest it's because the starting values aren't good, but I've tried tweaking them and I'm not sure that's the issue here. One post suggested using the function nlsLM from the minpack.lm package. Honestly I don't really know what that does, but I'll give it a try


library(minpack.lm)

nlsLM(Numeaten ~ I(k / 1 + exp(-r * (Temp - a))) + (1|Trial), data = Indata,
      start = list(k = log(50), r = 0.03, a = -5.7))

#Not working either. Volker suggested using SSasymp to generate starting parameters. His example code is in the next chunk. Giving it a try here

nls(Numeaten ~ SSasymp(Temp, Asym, R0, lrc), Indata)



#That's not working
# I just keep getting singularity errors. Here's some sample code Volker sent

####get porpostional survival
trisTS$propsurv<-trisTS$surv/100

#use simple nls fitting to get some first estimates for parameters for later model fitting. I fit it to each treatment first just to get general idea of how differen they were. You can ignore all of that.

#SSasymp and other functions are build in functions and you just need to define parameters, you can learn  more in help file.

#In your case you'd replace "ArrivalS" with temperature gradient treatment

asSH<-nls(propsurv~SSasymp(ArrivalS,Asym,R0,lrc),trisTS, subset=Nutrient=="H")
asSM<-nls(propsurv~SSasymp(ArrivalS,Asym,R0,lrc),trisTS, subset=Nutrient=="M")
asSN<-nls(surv~SSasymp(ArrivalS,Asym,R0,lrc),trisTS, subset=Nutrient=="L" )
gompM<-nls(surv~SSgompertz(ArrivalS,Asym,b2,b3),trisTS, subset=Nutrient=="M")
log4N<-nls(propsurv~SSfpl(ArrivalS,A,B,xmid,scal),trisTS,subset=Nutrient=="M")
summary(asSH)
summary(asSM)
summary(asSN)
AIC(asSM,gompM,log4N)




#with saturating/asymptotic regression mode allowing asymptote and slope to vary with Nutrient treatment (note you would just make R0 a constant in your case since you don't have different #treatments.
asymST<-nlme(propsurv~SSasymp(ArrivalS,Asym,R0,lrc),
             fixed=list(Asym+R0~Nutrient,lrc~1),
             random=Asym~1|Block,
             start=list(fixed=c(Asym=c(.6,.60,.25),R0=c(0.17,0.017,-0.04), lrc=-2)),
             weights = varComb(varIdent(form = ~1 | Nutrient)),
             data=trisTS)
AICc(asymST)
newdatsat<-expand.grid(Nutrient=levels(trisTS$Nutrient),ArrivalS=0:20)
newdatsat$pred<-predict(asymST,newdatsat,level=0)
ggplot(newdatsat,aes(x=ArrivalS,y=pred, color=Nutrient))+geom_line()+theme_bw()+geom_point(data=trisTS,aes(x=ArrivalS,y=propsurv, color=Nutrient))
anova(asymST)

#so lets fit the simple logistic function and get some stats. Noet that in your case Asym and mid are single falues, not a list. THey are a list here bc I have different starting values for each nutrient #treatment and fit parameters for each treatment

logST<-nlme(propsurv~SSlogis(ArrivalS,Asym,mid,scal),
            fixed=list(Asym+mid~Nutrient,scal~1),
            random=Asym~1|Block,
            start=list(fixed=c(Asym=c(0.6,0.4,0.2),mid=c(6,7,12),scal=2)),
            weights = varComb(varIdent(form = ~1 | Nutrient)),
            data=trisTS)
AICc(logST)
newdatlogST<-expand.grid(Nutrient=levels(trisTS$Nutrient),ArrivalS=0:20)
newdatlogST$pred<-predict(logST,newdatlogST,level=0)
ggplot(newdatlogST,aes(x=ArrivalS,y=pred, color=Nutrient))+geom_line()+theme_bw()+geom_point(data=trisTS,aes(x=ArrivalS,y=propsurv, color=Nutrient))
anova(logST)
###saturating gompertz function
gompST<-nlme(propsurv~SSgompertz(ArrivalS,Asym,b2,b3),
             fixed=list(Asym+b2~Nutrient, b3~1),
             random=Asym~1|Block,
             start=list(fixed=c(Asym=c(.60,.40,.20),b2=c(20,10,5),b3=.5)),
             weights = varComb(varIdent(form = ~1 | Nutrient)),
             data=trisTS)
newdatgompST<-expand.grid(Nutrient=levels(trisTS$Nutrient),ArrivalS=0:20)
newdatgompST$pred<-predict(gompST,newdatgompST,level=0)
ggplot(newdatgompST,aes(x=ArrivalS,y=pred, color=Nutrient))+geom_line()+theme_bw()+geom_point(data=trisTS,aes(x=ArrivalS,y=propsurv, color=Nutrient))


#Getting the predicted values with SE is proving to be really difficult. Here I'm just testing different methods


library(merTools)
feed.data <- as.data.frame(feed.data)

Buen.glmer <- feed.data%>%
    drop_na(Numeaten)%>%
    filter(Pred == "Buenoa")%>%
    glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp, 2)  +  scale(Predmass) + (1|Trial) + (1|Bath), family = "binomial", data = ., control = glmerControl(optimizer = "bobyqa"))


#Creating a dataset of the Buenoa data so I don't need to keep filtering
Buendata <- feed.data%>%
  drop_na(Numeaten)%>%
  filter(Pred == "Buenoa")


#Testing out predictInteral function for generating and graphing predicted values

#I was having an issue where I couldn't get the function to ignore/hold constant predmass, and it was giving me a jagged curve. Ultimately the only thing I could do that worked is to manually set the predmass values to the average value across all individuals of the same species. That's below:

#Set predator masses to the mean value
Buendata$Predmass <- mean(Buendata$Predmass, na.rm = T)

#Use predictInterval to calculate predictions and CI's then graph
Buendata%>%
  cbind(predictInterval(Buen.glmer, newdata = Buendata, type = "probability", which = "fixed",
                        include.resid.var = F))%>%  
  ggplot(aes(x = Temp, y = fit*100)) +
  geom_line() + 
  geom_point(aes(y = Numeaten)) +
  geom_ribbon(aes(ymax = upr*100, ymin = lwr*100), alpha = 0.5)


#Works!! Hot diggity. Moving on to the next chunk.




###Trying the bootMer function
#This function is just from the vignette for this function
mySumm <- function(.) { s <- sigma(.)
c(beta =getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta"))) }

#Here's the actual bootMer function
feed.data%>%
  na.omit()%>%
  filter(Pred == "Buenoa")%>%
  bootMer(Buen.glmer, mySumm, re.form = NA)






####Trying the marginaleffects package
library(marginaleffects)



#Predictions function to generate predictions. 
predictions(Buen.glmer, type = "response", variables = "Temp")
ggplot(aes(x = Temp))

#Trying the built in plotting function
plot_cap(Buen.glmer, condition= "Temp", type = "response")







###Trying ggeffects package
library(ggeffects)
ggpredict(Buen.glmer, terms = "Temp [all]")



#Here's the code for getting predicted values based on Volker's method in the heat tolerance data


#Create separate df of predicted values for each species. There's gotta be a more tidy way to do this but I can't think of one at the moment, barring a for loop

#Buenoa
Buenoa.predict <- feed.data%>%
  na.omit()%>%
  filter(Pred == "Buenoa")%>%
  mutate(Predicted = predictInterval(Buen.glmer, newdata = .)$fit,
         UprSE = predictInterval(Buen.glmer, newdata = .)$upr,
         LwrSE = Predicted - predictInterval(Buen.glmer, newdata = .)$lwr)

#Make new data frame to add each species' data to for graphing
Feed.fit.df <- data.frame(Trial = NA, Bath = NA, Location = NA,  Pred = NA, Temp = NA, Numeaten = NA, Predmass = NA, fit = NA, fit_link = NA, se_link = NA, fit_resp = NA, upr = NA, lwr = NA)


#make new data set for predictions
library(merTools)

Buen.predict <- feed.data%>%
  na.omit()%>%
  filter(Pred == "Buenoa")%>%
  predictInterval(Buen.glmer, newdata = ., which = "fixed", type = "probability")

#ndata <- feed.data%>%
#na.omit()%>%
#filter(Pred == "Buenoa")%>%
#with(data = ., data_frame(Temp = seq(min(Temp), max(Temp),
#length = 1000)))

#ndata<-add_column(ndata,fit=predictInterval(Buen.glmer, newdata=ndata, 
#which = "fixed"))

#get invers link function for prediction
ilink<-Buen.glmer@resp$family$linkinv
###add fit and se.fit on the link scale
Buen.predict<-bind_cols(Buen.predict, setNames(as_tibble(predictInterval(
  Buen.glmer, newdata = Buen.predict)[1:2]),
  c('fit_link','se_link')))
#create tge ubtercak and back transform
Indicaheat<-mutate(ndata,
                   fit_resp=ilink(fit_link),
                   upr=ilink(fit_link+(se_link)),
                   lwr=ilink(fit_link-(se_link)), 
                   Species = "Indica",
                   TempC = (Temp2 - 32)*(5/9))


