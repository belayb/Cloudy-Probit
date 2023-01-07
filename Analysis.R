# Load packages 
library(brms)
library(data.table)
library(dplyr)
library(tidybayes)
library(patchwork)
library(here)

#---------------Load data 
Pain_Data<-fread(here("Data/20190719_cloudy-data_motif-weather-bsl_794007-rows_10584_userids.csv"),header = T,sep=",")
length(unique(Pain_Data$UserId))# 10584 individuals and 794007 days of data read 
# Prepare Data for analysis 
# Focus- Impact of weather components (temprature, pressure, relative humidity and windspead) and pain, activity, mood and time spent outside
Pain_Data<-Pain_Data[!is.na(Pain_Data$painSeverity),]#452157 recoreds remain
length(unique(Pain_Data$UserId))# 10430 individuals remain
setDT(Pain_Data)

Pain_Data<-Pain_Data %>% 
  group_by(UserId) %>%
  dplyr::mutate(
    nentry = dplyr::n_distinct(day))
summary(Pain_Data$nentry)

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.0    74.0   156.0   181.5   286.0   449.0 

#---------------- We need at least two measurement per subject 
Pain_Data<-Pain_Data[Pain_Data$nentry>1,]#450311 records remain
length(unique(Pain_Data$UserId))# 8584 individuals remain

#---------------------- change calander time to process time -----------------------------------------------
Pain_Data <- Pain_Data %>%
  group_by(UserId) %>%
  mutate(Date=as.Date(day,format="%Y-%m-%d"),
         prev.entry_date = c(0, diff(Date)),
         time = cumsum(prev.entry_date))

length(unique(Pain_Data$UserId))
table(is.na(Pain_Data$meanTemp))#5.4% do not have weather info 

Pain_Data<-Pain_Data %>% 
  group_by(UserId) %>%
  dplyr::mutate(F_date = first(Date),L_date=last(Date),
                s_day = L_Date - F_date, 
                prop_entry = nentry/ as.numeric(s_day))%>%
  dplyr::filter(nentry>2)


#start model fitting 
Pain_Data$painSeverity<-ordered(Pain_Data$painSeverity)
Pain_Data$Mood1<-ifelse(Pain_Data$Mood>3,1,0)
Pain_Data$exercise1<-ifelse(Pain_Data$exercise>2,1,0)
Pain_Data$timeSpentOutside1<-ifelse(Pain_Data$timeSpentOutside>2,1,0)
Pain_Data$Belief<-ifelse(Pain_Data$belief>7,1,0)

Pain_Data$Belief<-factor(Pain_Data$Belief)
Pain_Data$Mood1<-factor(Pain_Data$Mood1)
Pain_Data$exercise1<-factor(Pain_Data$exercise1)
Pain_Data$timeSpentOutside1<-factor(Pain_Data$timeSpentOutside1)
Pain_Data<-Pain_Data[,c("UserId","Date", "time","Mood","painSeverity","exercise","timeSpentOutside","meanWindspeed","meanTemp",
                        "meanPressure", "meanRHX","belief","Belief","Mood1","exercise1","timeSpentOutside1",
                        "nentry")]

Base_dat<-fread(here("Data/20190703_10584_baselinewithoutcluster.csv"),header=T)#togetage

Base_dat<-Base_dat[,c("userid","sex","age")]
setDT(Base_dat)
setDT(Pain_Data)
names(Base_dat)[1]<-"UserId"

setkey(Base_dat,"UserId")
setkey(Pain_Data,"UserId")
Pain_Data<-Base_dat[Pain_Data]
# resclae humidity and pressure 
Pain_Data$meanRHX<-(Pain_Data$meanRHX)/10
Pain_Data$meanPressure<-(Pain_Data$meanPressure)/10

 # Fit model in brms

Model7<-brm(painSeverity ~ age+sex +Belief+ s(time,k=4)+meanTemp+meanPressure+meanRHX+meanWindspeed+Mood1+exercise1+ (1+meanTemp+meanPressure+meanRHX+meanWindspeed|UserId), data = Pain_Data,family=cumulative("probit"),
            set_prior("normal(0,1)",class="b"),set_prior("normal(0,1)",class="Intercept"),
            inits=0,chains = 2,core=2, iter = 6000,autocor=NULL)

# draw samples from the model 

p1<-Model7 %>%
  spread_draws(b_meanTemp, r_UserId[UserId,meanTemp]) %>%
  median_qi(UserId_mean = b_meanTemp + r_UserId)

p2<-Model7 %>%
  spread_draws(b_meanPressure, r_UserId[UserId,meanPressure]) %>%
  median_qi(UserId_mean = b_meanPressure + r_UserId)

p3<-Model7 %>%
  spread_draws(b_meanRHX, r_UserId[UserId,meanRHX]) %>%
  median_qi(UserId_mean = b_meanRHX + r_UserId)

p4<-Model7 %>%
  spread_draws(b_meanWindspeed, r_UserId[UserId,meanRHX]) %>%
  median_qi(UserId_mean = b_meanWindspeed + r_UserId)

# indentify non-overlaping CIs

p1$sig<-p1$.lower>0|p1$.upper<0
pd <- position_dodge(width=2)
#p1$UserId_mean_new<-ifelse(p1$UserId_mean>0, p1$.lower, p1$.upper)
p1<-left_join(p1, Pain_Data[!duplicated(Pain_Data[,"UserId"]),c("UserId","nentry")],by="UserId")
p1$length_95<-p1$.upper-p1$.lower
p1 <- p1 %>% mutate(ntimes=cut(nentry, breaks=c(-Inf,2, 7, 16, 40, 109, Inf),labels = c("2","3-7","8-16","17-40","41-109","110+")))

p2$sig<-p2$.lower>0|p2$.upper<0
p2<-left_join(p2, Pain_Data[!duplicated(Pain_Data[,"UserId"]),c("UserId","nentry")],by="UserId")
p2$length_95<-p2$.upper-p1$.lower
p2 <- p2 %>% mutate(ntimes=cut(nentry, breaks=c(-Inf,2, 7, 16, 40, 109, Inf),labels = c("2","3-7","8-16","17-40","41-109","110+")))

p3$sig<-p3$.lower>0|p3$.upper<0 
p3<-left_join(p3, Pain_Data[!duplicated(Pain_Data[,"UserId"]),c("UserId","nentry")],by="UserId")
p3$length_95<-p3$.upper-p1$.lower
p3 <- p3 %>% mutate(ntimes=cut(nentry, breaks=c(-Inf,2, 7, 16, 40, 109, Inf),labels = c("2","3-7","8-16","17-40","41-109","110+")))

p4$sig<-p4$.lower>0|p4$.upper<0 
p4<-left_join(p4, Pain_Data[!duplicated(Pain_Data[,"UserId"]),c("UserId","nentry")],by="UserId")
p4$length_95<-p4$.upper-p1$.lower
p4 <- p4 %>% mutate(ntimes=cut(nentry, breaks=c(-Inf,2, 7, 16, 40, 109, Inf),labels = c("2","3-7","8-16","17-40","41-109","110+")))


p11<-p1%>%
  ggplot(aes(y = reorder(UserId,UserId_mean), x = UserId_mean, xmin = .lower, xmax = .upper,color=sig)) +
  geom_pointintervalh(position = pd, size=1)+
  labs(title = "Temperature")+
  theme(panel.grid.minor = element_blank(),legend.position = "none",
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_color_manual(values=c("gray","blue")) + 
  geom_vline(xintercept = 0) +
  xlab("95% CI")+
  ylab("Participants")+xlim(-0.2,0.2)+
  geom_vline(xintercept = -0.003, linetype='dashed', color='red') +
  coord_flip()

p12<-p1 %>%  
  group_by(ntimes,sig) %>% 
  summarise(CIlength = mean(length_95)) %>% 
  ggplot(aes(ntimes, CIlength, colour=sig, group=sig)) + 
  geom_line(size=1.5)+ylim(0,0.25)+
  scale_color_manual(values=c("#A9A9A9","blue"))+
  xlab("Number of repeated observation")+ylab("95% CI length")+
  theme_bw()+ theme(legend.position = "none")

p21<-p2%>%
  ggplot(aes(y = reorder(UserId,UserId_mean), x = UserId_mean, xmin = .lower, xmax = .upper,color=sig)) +
  geom_pointintervalh(position = pd, size=1)+
  labs(title = "Relative humidity")+
  theme(panel.grid.minor = element_blank(),legend.position = "none",
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_color_manual(values=c("gray","blue")) + 
  geom_vline(xintercept = 0) +
  xlab("95% CI")+
  ylab("Participants")+
  geom_vline(xintercept = 0.041, linetype='dashed', color='red') + xlim(-0.5, 0.5)+
  coord_flip()

p22<-p2 %>%  
  group_by(ntimes,sig) %>% 
  summarise(CIlength = mean(length_95)) %>% 
  ggplot(aes(ntimes, CIlength, colour=sig, group=sig)) + 
  geom_line(size=1.5)+ylim(0,0.7)+
  scale_color_manual(values=c("#A9A9A9","blue"))+
  xlab("Number of repeated observation")+ylab("95% CI length")+
  theme_bw()+ theme(legend.position = "none")


p31<-p3%>%
  ggplot(aes(y = reorder(UserId,UserId_mean), x = UserId_mean, xmin = .lower, xmax = .upper,color=sig)) +
  geom_pointintervalh(position = pd, size=1)+
  labs(title = "Pressure")+
  theme(panel.grid.minor = element_blank(),legend.position = "none",
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_color_manual(values=c("gray","blue")) + 
  geom_vline(xintercept = 0) +
  xlab("95% CI")+
  ylab("Participants")+xlim(-0.4,0.4)+
  geom_vline(xintercept = -0.01, linetype='dashed', color='red') +
  coord_flip()

p32<-p3 %>%  
  group_by(ntimes,sig) %>% 
  summarise(CIlength = mean(length_95)) %>% 
  ggplot(aes(ntimes, CIlength, colour=sig, group=sig)) + 
  geom_line(size=1.5)+ylim(0,0.4)+
  scale_color_manual(values=c("#A9A9A9","blue"))+
  xlab("Number of repeated observation")+ylab("95% CI length")+
  theme_bw()+ theme(legend.position = "none")


p41<-p4%>%
  ggplot(aes(y = reorder(UserId,UserId_mean), x = UserId_mean, xmin = .lower, xmax = .upper,color=sig)) +
  geom_pointintervalh(position = pd, size=1)+
  labs(title = "Wind speed")+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_color_manual(values=c("gray","blue")) + 
  geom_vline(xintercept = 0) +
  xlab("95% CI")+
  ylab("Participants")+xlim(-0.2,0.2)+
  geom_vline(xintercept = 0.012, linetype='dashed', color='red') +
  coord_flip()
p42<-p4 %>%  
  group_by(ntimes,sig) %>% 
  summarise(CIlength = mean(length_95)) %>% 
  ggplot(aes(ntimes, CIlength, colour=sig, group=sig)) + 
  geom_line(size=1.5)+ylim(0,0.4)+
  scale_color_manual(values=c("#A9A9A9","blue"))+
  xlab("Number of repeated observation")+ylab("95% CI length")+
  theme_bw()+ theme(legend.position = "none")

#-------------Exposure hetroginity plot main paper 

p_all1<-p11+p12+p21+p22+p31+p32+p41+p42 + plot_layout(widths = c(1.5, 1))


# Do conditions expalain class membership - outcome (low senstive, high senstative, undetrmined)
cond_data<-Base_dat[Base_dat$userid%in%Pain_Data[!duplicated(Pain_Data[,"UserId"]),]$UserId,
                    c(4:13, 21)]

cond_data$number_cond<-rowSums(cond_data[, 2:10])
cond_data_long <- gather(cond_data, condition, presence, condsum.ra:cond.final.other, factor_key=TRUE)
cond_data_long<-cond_data_long[!(cond_data_long$presence==0),]
cond_data_long<-cond_data_long %>% 
  mutate(Chronic_pain_type = case_when(
    condition%in%c("condsum.ra","condsum.gou","condsum.spa") ~ 'Inflammatory arthritis pain',
    condition%in%c("condsum.oa","condsum.artr.unsp","condsum.cwpfm") ~ 'Non-Inflammatory MSK pain',
    TRUE ~ 'Other Chronic pain ') )

cond_data_long<-cond_data_long%>%
  group_by(userid)%>%
  mutate(n.condition=sum(presence))
cond_data_long<-cond_data_long%>%
  group_by(userid)%>%
  mutate(n.pain_type=n_distinct(Chronic_pain_type))

with(cond_data_long, table(n.pain_type))
length(unique(cond_data_long[cond_data_long$n.pain_type==1,]$userid))
cond_data_long<-cond_data_long[cond_data_long$n.pain_type==1,c("userid","Chronic_pain_type")]
cond_data_long<-cond_data_long[!duplicated(cond_data_long[,c("userid","Chronic_pain_type")]),]

length(unique(cond_data_long$userid))
names(cond_data_long)[1]<-"UserId"
# windspeed effect
cond_data_long_windspeed<-left_join(cond_data_long, p4, by="UserId")
head(cond_data_long_windspeed)
#Define the multinomial outcome
cond_data_long_windspeed<-cond_data_long_windspeed %>% 
  mutate(Outcome_Class_type = case_when(
    .lower>0 ~ 'High value sensitive',
    .upper<0 ~ 'Low value sensitive',
    TRUE ~ 'Undetermined') )

with(cond_data_long_windspeed, table(Chronic_pain_type,Outcome_Class_type) )
colors<-c("#a6a6a6","#E9C46A","#2F575D")#,  "#f67e7d", "#843b62"
bar_windspeed<-cond_data_long_windspeed%>%
  group_by(Outcome_Class_type,Chronic_pain_type)%>%
  summarise(n = n()) %>%
  mutate(Prevalence = 100*(n / sum(n)))%>%
  #filter(!(Outcome_Class_type=="Undetermined"))%>%
  ggplot(aes(x=Outcome_Class_type, y=Prevalence, fill=forcats::fct_rev(Chronic_pain_type))) +
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values=colors)+
  labs(fill = "Chronic pain class",title = "Windspeed",element_text(color = "blue", size = 10))+
  guides(guide_legend(title.position="top", title.hjust = 0.5))+
  ylab("Percent occurrence")+
  xlab("Cluster")+
  theme_bw()

# Temperature effect
cond_data_long_temperature<-left_join(cond_data_long, p1, by="UserId")
head(cond_data_long_temperature)
#Define the multinomial outcome
cond_data_long_temperature<-cond_data_long_temperature %>% 
  mutate(Outcome_Class_type = case_when(
    .lower>0 ~ 'High value sensitive',
    .upper<0 ~ 'Low value sensitive',
    TRUE ~ 'Undetermined') )
with(cond_data_long_temperature, table(Chronic_pain_type,Outcome_Class_type) )
bar_temp<-cond_data_long_temperature%>%
  group_by(Outcome_Class_type,Chronic_pain_type)%>%
  summarise(n = n()) %>%
  mutate(Prevalence = 100*(n / sum(n)))%>%
  #filter(!(Outcome_Class_type=="Undetermined"))%>%
  ggplot(aes(x=Outcome_Class_type, y=Prevalence, fill=forcats::fct_rev(Chronic_pain_type))) +
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values=colors)+
  labs(fill = "Chronic pain class", title="Temperature",element_text(color = "blue", size = 10))+
  guides(guide_legend(title.position="top", title.hjust = 0.5))+
  ylab("Percent occurrence")+
  xlab("Cluster")+
  theme_bw()
# relative humidity effect
cond_data_long_RH<-left_join(cond_data_long, p2, by="UserId")
head(cond_data_long_RH)
#Define the multinomial outcome
cond_data_long_RH<-cond_data_long_RH %>% 
  mutate(Outcome_Class_type = case_when(
    .lower>0 ~ 'High value sensitive',
    .upper<0 ~ 'Low value sensitive',
    TRUE ~ 'Undetermined') )
with(cond_data_long_RH, table(Chronic_pain_type,Outcome_Class_type) )
bar_RH<-cond_data_long_RH%>%
  group_by(Outcome_Class_type,Chronic_pain_type)%>%
  summarise(n = n()) %>%
  mutate(Prevalence = 100*(n / sum(n)))%>%
  #filter(!(Outcome_Class_type=="Undetermined"))%>%
  ggplot(aes(x=Outcome_Class_type, y=Prevalence, fill=forcats::fct_rev(Chronic_pain_type))) +
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values=colors)+
  labs(fill = "Chronic pain class", title="Relative humidity", element_text(color = "blue", size = 10))+
  guides(guide_legend(title.position="top", title.hjust = 0.5))+
  ylab("Percent occurrence")+
  xlab("Cluster")+
  theme_bw()

# Pressure effect
cond_data_long_pressure<-left_join(cond_data_long, p3, by="UserId")
head(cond_data_long_pressure)
#Define the multinomial outcome
cond_data_long_pressure<-cond_data_long_pressure %>% 
  mutate(Outcome_Class_type = case_when(
    .lower>0 ~ 'High value sensitive',
    .upper<0 ~ 'Low value sensitive',
    TRUE ~ 'Undetermined') )
with(cond_data_long_pressure, table(Chronic_pain_type,Outcome_Class_type) )
bar_pressure<-cond_data_long_pressure%>%
  group_by(Outcome_Class_type,Chronic_pain_type)%>%
  summarise(n = n()) %>%
  mutate(Prevalence = 100*(n / sum(n)))%>%
  #filter(!(Outcome_Class_type=="Undetermined"))%>%
  ggplot(aes(x=Outcome_Class_type, y=Prevalence, fill=forcats::fct_rev(Chronic_pain_type))) +
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values=colors)+
  theme(axis.text.x = element_text(color = "black",face = "plain"))+
  labs(fill = "Chronic pain class", title="Pressure",element_text(color = "blue", size = 10))+
  guides(guide_legend(title.position="top", title.hjust = 0.5))+
  ylab("Percent occurrence")+
  xlab("Cluster")+
  theme_bw()

(bar_temp + bar_RH) /(bar_pressure+bar_windspeed) +plot_layout(guides = "collect")+
  guides(guide_legend(title.position="top", title.hjust = 0.5))&theme(legend.position = "bottom",
                                                                      legend.title = element_text(face = "bold", size = 12))



# ------------------- Tabulation of diffrent wether sensetive individuals --------------



p1<-p1 %>% 
  mutate(Outcome_Class_type = case_when(
    .lower>0 ~ 'High value sensitive',
    .upper<0 ~ 'Low value sensitive',
    TRUE ~ 'Undetermined') )

p2<-p2 %>% 
  mutate(Outcome_Class_type = case_when(
    .lower>0 ~ 'High value sensitive',
    .upper<0 ~ 'Low value sensitive',
    TRUE ~ 'Undetermined') )

p3<-p3 %>% 
  mutate(Outcome_Class_type = case_when(
    .lower>0 ~ 'High value sensitive',
    .upper<0 ~ 'Low value sensitive',
    TRUE ~ 'Undetermined') )

p4<-p4 %>% 
  mutate(Outcome_Class_type = case_when(
    .lower>0 ~ 'High value sensitive',
    .upper<0 ~ 'Low value sensitive',
    TRUE ~ 'Undetermined') )

names(p1)[2]<-"Weather";names(p2)[2]<-"Weather";names(p3)[2]<-"Weather";names(p4)[2]<-"Weather"
p_all<-rbind(p1[,c("UserId","Weather","Outcome_Class_type")],p2[,c("UserId","Weather","Outcome_Class_type")],
             p3[,c("UserId","Weather","Outcome_Class_type")],p4[,c("UserId","Weather","Outcome_Class_type")])

head(p_all)
p_all_wide<-spread(p_all, Weather, Outcome_Class_type)
head(p_all_wide)
p_all_wide<-p_all_wide %>% 
  mutate(P = case_when(
    meanPressure=="Undetermined" ~ 'U',
    meanPressure=="Low value sensitive" ~ 'L',
    TRUE ~ 'H') )
p_all_wide<-p_all_wide %>% 
  mutate(R = case_when(
    meanRHX=="Undetermined" ~ 'U',
    meanRHX=="Low value sensitive" ~ 'L',
    TRUE ~ 'H') )
p_all_wide<-p_all_wide %>% 
  mutate(TT = case_when(
    meanTemp=="Undetermined" ~ 'U',
    meanTemp=="Low value sensitive" ~ 'L',
    TRUE ~ 'H') )
p_all_wide<-p_all_wide %>% 
  mutate(W = case_when(
    meanWindspeed=="Undetermined" ~ 'U',
    meanWindspeed=="Low value sensitive" ~ 'L',
    TRUE ~ 'H') )

p_all_wide$PRTW<-paste0(p_all_wide$P,p_all_wide$R, p_all_wide$TT,p_all_wide$W)

PRTW_count<-p_all_wide%>%
  group_by(PRTW)%>%
  summarise(n = n())

PRTW_plot<-p_all_wide%>%
  group_by(PRTW)%>%
  summarise(n = n()) %>%
  mutate(Prevalence = 100*(n / sum(n)))%>%
  filter(!(PRTW=="UUUU"))%>%
  ggplot(aes(x= reorder(PRTW,-Prevalence), y=Prevalence)) +
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values=colors)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Weather sensetivity combinations")