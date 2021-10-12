
rm(list = ls(all.names = TRUE))
gc()

library(tidyverse)
library(glmnet)
library(survival)
library(ggRandomForests)
library(randomForestSRC)
library(ggplot2)

library(tidymodels)
library(recipes)
library(modeldata)
library(themis)

library(gtsummary)
library(biostat3)
library(survminer)
library(broom)
library(cowplot)

path = c("//home.gu.gu.se/home-XB$/xbocda/Documents/RCVR/dat")
path <- c("<insert path here>")    # path
A <- readRDS(paste0(path, "<Name of data set>.rds"))    # Name of data set

#A$Adverse.Events <-  factor(A$Adverse.Events , labels = c("No", "Yes") )
#A$MKR <-  factor(A$MKR , labels = c("Radikalitet Nej", "Radikalitet Ja") )
#A$Op.typ  <-  factor(A$Op.typ  , labels = c("Resektion" , "APE/ELAPE", "Hartmann" ,  "Övrigt" ) )
#A$Recurrence  <-  factor(A$Recurrence  , labels = c("Recurrence: No" , "Recurrence:  Yes") )
#A$SKOLJ  <-  factor(A$SKOLJ  , labels = c("No" , "Yes" , "Not applicable") )
# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
# To increase speed, we can fit the models in parallel :
library("doFuture")
all_cores <- parallel::detectCores(logical = FALSE) - 1

registerDoFuture()
cl <- makeCluster(all_cores)
plan(future::cluster, workers = cl)

# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

A2 <- A %>% mutate(EVENT2 = case_when(  EVENT==0  ~ 1 ,
                                        EVENT==1  ~ 0 ,
                                        EVENT==2  ~ 2   ) )


#  0/1/2 for censored, recurrence, dead
#   0   1   2
#  747  78 383
A2$EVENT_  <- factor(A2$EVENT2 , labels = c("Censored" , "Recurrence" , "Deceased") )
A2$time <- A2$FU

# Out of 1208 patienter 1129 were either right censored at follow-up (746 st) or dead (383)
#EVENT    min       median   max
#Rec      0.110     1.71    8.21
#Cens     5.00      6.83   9.74
#Död      0.00548   3.34   9.34

A2 %>% group_by(EVENT_) %>%
  summarise(min = min(FU, na.rm=T), median = median(FU, na.rm=T), max = max(FU, na.rm=T))
#######################################################################################################################
library(prodlim)
library(cmprsk)
library(riskRegression)
library(Publish)
library(timereg)

#  0/1/2 for censored, recurrence, dead
#   0   1   2
#  746  79 383

#   Cause-specific Cox proportional hazard regression
B <- A2 %>% dplyr::select( c(FU , EVENT2 ,
                    Adverse.Events ,
                    MKR ,
                    N.PAD ,
                    Neoadjuv.beh ,
                    Op.typ ,
                    pT ,
                    SKOLJ ,
                    Perop.rektskölj ) )


B$Adverse.Events <-  factor(B$Adverse.Events , labels = c("No", "Yes") )
B$MKR <-  factor(B$MKR , labels = c("Radikalitet Nej", "Radikalitet Ja") )
B$Op.typ  <-  factor(B$Op.typ  , labels = c("Resektion" , "APE/ELAPE", "Hartmann"  ) )
B$SKOLJ  <-  factor(B$SKOLJ  , labels = c("No" , "Yes" , "Not applicable") )
B$N.PAD  <-  factor(B$N.PAD  , labels = c("No Nodes" , "Yes Nodes") )
B$pT  <-  factor(B$pT  , labels = c("Minor" , "Modest/Major") )
B$Neoadjuv.beh <-  factor(B$Neoadjuv.beh  , labels = c("No" , "Yes" ) )

################################################################################################################################
CX8 <- CSC(Hist(FU , EVENT2) ~  Adverse.Events + MKR + N.PAD + Neoadjuv.beh  + Op.typ + pT + SKOLJ , d = B , cause = 1)
FGR8 <- FGR(Hist(FU , EVENT2) ~  Adverse.Events + MKR + N.PAD + Neoadjuv.beh  + Op.typ + pT + SKOLJ  , d = B , cause = 1)
CX81 <- publish(CX7) ;
CX81$`Cause 1`
z<- summary(FGR8) ; z$coef ; z$conf.int ;

#######################################################################################################################
#  0/1/2 for censored, recurrence, dead
#   0   1   2
#  746  78 383
out1 <- comp.risk(Event(FU, EVENT2 ) ~  1, data = B, cause = 1)
summary(out1)
pout1 <- predict(out1, X = 1,  times = c(5, 10))
names(pout1)
pout1$P1 ;
pout1$se.P1 ;
LCL = pout1$P1 - 1.96* pout1$se.P1 ;
UCL = pout1$P1 + 1.96* pout1$se.P1 ;

glimpse(B)

B$EVENT_  <- factor(B$EVENT2 , labels = c("Censored" , "Recurrence" , "Deceased") )

# Cumulative incidence
ci_fit <-
  cuminc(
    ftime = B$FU,
    fstatus =A2$EVENT_,
    cencode = "Censored"  )

print(ci_fit)


plt_ <-  ggcompetingrisks(  ci_fit , conf.int = TRUE , ggtheme = theme_survminer())

plt_2 <- plt_    +
  theme(legend.position = c(0.175, 0.75)) +
  scale_linetype_manual(values=c("dotted"))+
  scale_color_manual(values=c("red", "black"))+
  scale_size_manual(values=c(1, 2.5))+
  geom_line(     size=1.0)   +
  theme(legend.title = element_text(colour="white")) +
  theme(legend.text = element_text(size = 16, colour = "black", angle = 0)) +
  theme(legend.key.size = unit(0.75, "cm")) +
  theme(panel.background = element_rect(fill = "white"))  +
  theme(plot.background = element_rect()) +
  xlab("Years") + ylab("Cause-specific cumulative incidence") +
  theme(plot.title = element_text(size = 0 , colour = "white")) +
  theme( strip.background = element_blank(),  strip.text.x = element_blank() ) +
  annotate("text",
           x = c(4.25, 7.9), y = c(0.115, 0.115), # x and y coordinates of the text
           label = c("5.8%(4.5; 7.1)" , "6.9%(5.3; 8.4)" ) ,
           size = 4.9, hjust = 0)

mel_fit <- survfit(
  Surv(FU, ifelse(EVENT2 != 2, 1, 0)) ~ 1,
  data = B)

num <- ggsurvplot(
  fit = mel_fit,
  risk.table = TRUE,
  risk.table.y.text = FALSE,
  ylab = "Years",
  risk.table.fontsize = 4.5,
  palette = c("white") ,
  tables.theme = theme_survminer(font.main = 12 ))

mel_plot <- plt_2

COWW <-  cowplot::plot_grid(
  mel_plot,
  num$table + theme_cleantable(),
  nrow = 2,
  rel_heights = c(4, 1),
  align = "v",
  axis = "b" )

dev.new(width=20, height=20)
window()
pd <- position_dodge(.25)
par(mar = c(0.4, 0.1, 0.4, 0.1))

COWW



#######################################################################################################################
pout1 <- predict(out1, X = 1,  times = c(5, 10))


out2 <- comp.risk(Event(FU, EVENT2 ) ~  Adverse.Events, data = B, cause = 1)
summary(out2)
names(out2)
ndata<-data.frame(Adverse.Events=c(0, 1))
pout2 <- predict(out2, newdata = ndata ,   times = c( 5, 10))

pout2$P1 ;
LCL = pout2$P1 - 1.96* pout2$se.P1 ;
UCL = pout2$P1 + 1.96* pout2$se.P1 ;
LCL ; UCL ;


#  0/1/2 for censored, recurrence, dead

# Cumulative incidence
ci_ <-
  cuminc(
    ftime = A2$time,
    fstatus =A2$EVENT2,
    group = A2$Adverse.Events ,
    cencode = 0  )

ciplotdat1 <-
  ci_ %>%
  list_modify("Tests" = NULL) %>%
  map_df(`[`, c("time", "est"), .id = "id") %>%
  filter(id %in% c("0 1", "1 1")) %>%
  mutate(Adverse.Events = recode(
    id,
    "0 1" = "No adverse event",
    "1 1" = "Adverse event")
  )




#   Get a plot from base R, ggcompetingrisks, or ggplot
mel_plot <-
  ggplot(ciplotdat1, aes(x = time, y = est, color =  Adverse.Events )) +
  geom_step(lwd = 1.2)  +
  ylim(c(0, 0.3)) +                             # 0.5
  coord_cartesian(xlim = c(0, 10)) +
#  scale_x_continuous(breaks = seq(0, 2, 4, 6, 8)) +
  theme_classic() +
  theme(plot.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom") +

                               theme(legend.text = element_text(size = 12, colour = "black", angle = 0)) +
                                theme(legend.position = c(0.175, 0.4)) +

                                 theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12)) +

  labs(x = "Years",
       y = "Cumulative incidence") +
  annotate("text", x = 0.075, y = 0.25, hjust = 0, size = 4.5,
           label = paste0(
             "p-value : ",
             ifelse(ci_$Tests[1, 2] < .001,
                    "<.001",
                    round(ci_$Tests[1, 2], 3)))) +
  annotate("text",
           x = c(4.25, 7.9 , 4.25, 7.9), y = c(0.115, 0.115 , 0.0475, 0.055), # x and y coordinates of the text
           label = c("8.6%(6.4; 10.7)" , "9.6%(7.3; 11.8)" , "2.4%(1.1; 3.7)" , "3.6%(1.6; 5.6)" ) ,
           size = 4.5, hjust = 0)

#   Get the number at risk table from a ggsurvplot using the survfit where all events count as a single composite endpoint

                              ##########
                              TMP <- A2 %>% mutate( Adverse.Events_ = if_else(Adverse.Events==1, 0, 1) )

mel_fit <- survfit(
  Surv(FU, ifelse(EVENT2 != 2, 1, 0)) ~ Adverse.Events_,
  data = TMP )

num <- ggsurvplot(     fit = mel_fit,
  risk.table = TRUE,
  risk.table.y.text = FALSE,
  ylab = "Years",
  risk.table.fontsize = 4.2,
  tables.theme = theme_survminer(font.main = 10),
  title = "Test"
)

COWW2 <-  cowplot::plot_grid(
  mel_plot,
  num$table + theme_cleantable(),
  nrow = 2,
  rel_heights = c(4, 1),
  align = "v",
  axis = "b" )

COWW2
