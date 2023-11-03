library(tidyverse)
#Import data
datLog <- read.csv("Data/Data.logger2.csv")

#Reformat the date time variable
datLog <- mutate(datLog, date.time = mdy_hms(datLog$date.time))



#Calculate daily min, max, and avg temps
finalDf <- datLog%>%
  filter(X. %in% 1:4000 | X. %in% 7500:20150 | X. %in% 22500:24400 | X. %in% 25400:31760 | X. > 33800)%>%#remove dates logger wasn't submerged
  group_by(as.Date(date.time))%>%
  summarize(minTemp = min(temp),
            maxTemp = max(temp),
            meanTemp = mean(temp))%>%
  rename(date = "as.Date(date.time)")


#Plot
finalDf%>%
  #filter(X. %in% 1:4000 | X. %in% 7500:20150 | X. %in% 22500:24400 | X. %in% 25400:31760 | X. > 33800)%>%
  #filter(X. %in% 24400:25400) %>%
  ggplot(aes(x = date,y = meanTemp)) +
  #geom_point(size = 0.04) +
  # geom_smooth(color = "#FFC107", fill = "#FFC107", alpha = 0.2) +
  # geom_smooth(aes(y = minTemp), color = "#1E88E5", fill = "#1E88E5", alpha = 0.2) +
  # geom_smooth(aes(y = maxTemp), color = "#D81B60", fill = "#D81B60", alpha = 0.2) +
  geom_smooth(color = "black", se = F) +
  geom_smooth(aes(y = minTemp), color = "black", linetype = "dashed", se = F) +
  geom_smooth(aes(y = maxTemp), color = "black", linetype = "dashed", se = F) +
  theme_classic() +
  labs(x = "Date", y = "Temperature (°C)") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        plot.margin = margin(t = 10, b = 5, l = 10, r = 25)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b")

ggsave("Figures/dataLogger.jpeg", width = 6.5, height = 5.6)

 