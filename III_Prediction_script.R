
gc()
options("install.lock"=FALSE)
options(scipen = 6, digits = 4) # view outputs in non-scientific notation

library(tidyverse)
library(glmnet)
library(glmnetUtils)
library(survival)
library(ggRandomForests)
library(randomForestSRC)
library(ggplot2)
library(Hmisc)
library(rms)
library(tidymodels)
library(recipes)
library(modeldata)
library(themis)
library(purrr)
library(gtsummary)
library(biostat3)
library(survminer)
library(broom)
library(cowplot)
library(rsample)

# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
# To increase speed, we can fit the models in parallel :
library("doFuture")
library(parallel)
library(doParallel)
all_cores <- parallel::detectCores(logical = FALSE) - 1

registerDoFuture()
cl <- makeCluster(all_cores/2)
plan(future::cluster, workers = cl)
registerDoParallel(cluster)


# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#
path <- c("<insert path here>")    # Path
A <- readRDS(paste0(path, "<Name of data set>.rds"))   # Name of data set
B <- A %>% mutate(Recurrence = case_when(  EVENT2==1 ~ 1 , EVENT2==0 ~ 0 , EVENT2==2 ~ 0 ) ,
                  Recurrence = factor(Recurrence)) %>%
  dplyr::select(FU, Recurrence , Adverse.Events , MKR , N.PAD , Op.typ ,
                Neoadjuv.beh , Op.typ , SKOLJ , pT , Perop.rektskölj )


B <- B %>% mutate(Op.typ2 = case_when( Op.typ=="Resektion" &  Perop.rektskölj==0 ~ 1 ,   
                                       Op.typ=="Resektion" & Perop.rektskölj==1 ~ 2 ,  
                                       Op.typ=="Hartmann" &  Perop.rektskölj==0 ~ 3 , 
                                       Op.typ=="Hartmann" & Perop.rektskölj==1 ~ 4 ,  
                                       Op.typ=="APE/ELAPE"  ~ 5 ))    %>%
  mutate( Op.typ_ = factor( Op.typ2 ),
          Op.typx = factor(Op.typ2, levels = c(1, 2, 3, 4, 5),
                           labels = c("Resection_w.o_wash",
                                      "Resection_w_wash",
                                      "Hartmann_w.o_wash",
                                      "Hartmann_w_wash",
                                      "APE.ELAPE")))

REC <- recipe( Recurrence ~ . , data = B) %>%
  step_impute_knn( Adverse.Events , MKR , N.PAD , Neoadjuv.beh , Op.typ , Op.typ_, Op.typx, pT ) %>%
  step_dummy( Adverse.Events , MKR , N.PAD ,  Neoadjuv.beh  , pT )

prepped_recipe <- prep(REC )
Bprep  <- juice(prepped_recipe)
names(Bprep)

C <- REC %>%
  # apply the recipe to the training data
  prep(B) %>%
  # extract the pre-processed training dataset
  juice()
#_____________________________________________________________________________________
C <- C %>% rename(
  Adverse.Events_Yes = Adverse.Events_Yes ,
  Radical.surgery_No.vs.Yes = MKR_Radikalitet.Ja,
  Neoadj.treatment_Yes.vs.No = Neoadjuv.beh_Yes ,
  Nodes_Yes.vs.No = N.PAD_Yes.Nodes  ,
  Tumor.stage_3to4.vs.0to2 = pT_Modest.Major ) %>%
  dplyr::select(  c(FU, Recurrence ,
                    Adverse.Events_Yes ,
                    Radical.surgery_No.vs.Yes ,
                    Neoadj.treatment_Yes.vs.No ,
                    Op.typ_ ,
                    Op.typx,
                    Nodes_Yes.vs.No ,
                    Tumor.stage_3to4.vs.0to2 )) %>%
  mutate(EVENT = case_when( Recurrence == 0 ~ 0 , Recurrence == 1  ~ 1   )) %>%
  mutate(EVENT = factor(EVENT)) %>% dplyr::select(-c(Recurrence))

C2 <- tibble(time=survival::Surv(C$FU, C$EVENT, typ = "right")[,1],
             status=survival::Surv(C$FU, C$EVENT, typ = "right")[,2]) %>%
  cbind(C)
C2_0 <- C2 %>% mutate(TRT = Op.typx)
#_____________________________________________________________________________________
# Internal validation. Results  
library(hdnom)
time <- as.numeric(C2_0$FU)
event <- as.numeric(C2_0$status)
Y_ <- survival::Surv(time, event, typ = "right")
X_ <- model.frame(TRT ~ Adverse.Events_Yes +
                    Radical.surgery_No.vs.Yes +
                    Neoadj.treatment_Yes.vs.No +
                    Nodes_Yes.vs.No +
                    Tumor.stage_3to4.vs.0to2 +
                    TRT, data = C2_0)
X_ <- glmnet::makeX(X_)
#_____________________________________________________________________________________
fit <- fit_lasso(X_, Y_, nfolds = 10, rule = "lambda.1se", seed = 123)
lambda <- fit$lambda
nom <- as_nomogram(
  fit, X_, time, event,
  pred.at = 5)

val_int <- validate( X_, time, event,
                     model.type = "lasso",
                     lambda = lambda, alpha = 1,
                     method = "bootstrap", boot.times = 1000,
                     tauc.type = "UNO", tauc.time = 5, 
                     seed = 58902, trace = FALSE)
val_int2 <- summary(val_int)
#_____________________________________________________________________________________

cal_int <- calibrate( X_, time, event, model.type = "lasso",
                      lambda = lambda, alpha = 1,
                      method = "bootstrap", boot.times = 1000,
                      pred.at = 8, ngroup = 3,  seed = 412, trace = FALSE )

cl<-summary(cal_int)
Cl <- tibble( Pr=1-cl[,1], 
              Ob=1-cl[,2],
              U=1-cl[,3], 
              L=1-cl[,4])

fig2 <-    ggplot(Cl, aes(x=Pr, y=Ob)) +
  geom_line(colour = "black" , size = 0.5, lty = 1) + 
  geom_point(colour ="black", size = 2,  ) + 
  geom_errorbar(aes(ymin=L, ymax=U),colour = "black",  width=.005, size=0.5, lty=1)+
  scale_x_continuous(limits=c(0, 0.225))+
  scale_y_continuous(limits=c(0, 0.225))+
  geom_abline(intercept=0, slope=1, color='grey', size=0.5, lty=1)+
  theme_bw() +
  xlab("Observed cumulative five year recurrence risk") +
  ylab("Predicted cumulative five year recurrence risk") +  
  theme(axis.title.x = element_text(color="black", size=12),   
        axis.title.y = element_text(color="black", size=12),   
        axis.text.x = element_text(color = "black", size = 11,
                                   angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.text.y = element_text(color = "black", size = 11,
                                   angle = 0, hjust = 1, vjust = 0, face = "plain"))  

# The median (Q1; Q3) time-dependent AUC at five year follow-up is 0.85(0.84; 0.86). 

library(rms)
library(Hmisc)

pr <- cph(Y_ ~ Radical.surgery_No.vs.Yes +
            Neoadj.treatment_Yes.vs.No +
            Nodes_Yes.vs.No +
            Tumor.stage_3to4.vs.0to2 +
            TRT, data = C2_0, x=TRUE, y=TRUE)
po <- cph(Y_ ~  Adverse.Events_Yes + 
            Radical.surgery_No.vs.Yes +
            Neoadj.treatment_Yes.vs.No +
            Nodes_Yes.vs.No +
            Tumor.stage_3to4.vs.0to2 +
            TRT, data = C2_0, x=TRUE, y=TRUE)


lrtest(po, pr)   #   Indicates that Adverse events adds value. To what extent?

pr   #R2 = 0.162
po   #R2 = 0.179
cREDUCED <- ((pr$stats['Dxy'])/2)+0.5
cFULL  <- ((po$stats['Dxy'])/2)+0.5

lrREDUCED <- pr$stats['Model L.R.']
lrFULL <- po$stats['Model L.R.']
rREDUCED  <- pr$stats['R2']
rFULL  <- po$stats['R2']

#LR Chi2 Pre =  120.86 
#LR Chi2 Post = 134.46 
# Adequacy of base model
lrFULL/ lrREDUCED   # 0.8989
# Fraction of new information from AE:
#  fraction of diagnostic information that is new information from AE
lrFULL/ lrREDUCED -1  # 11%

#  Fraction of new information” is the proportion of total predictive information of the model
# including AE that was added by AE.
# one minus the ratio of the variance of pre-test probability of disease  to the variance of post-test
#   More discriminating models provide a greater variety of predictions, subject to assuming the model is well-calibrated
#   AE adds a fraction of 0.10 of new diagnostic information to the other variables.

fig22 <- fig2 +  annotate("text",
                  y = c(0.08 , 0.08, 0.05, 0.05, 0.021),    #  c(0.22 , 0.22, 0.19, 0.19, 0.16)
                  x =  c(0.12 , 0.185, 0.125, 0.185 , 0.12), 
                  label =  c("Adverse events\nnot included:" ,
                             "Adverse events \nincluded:",
                             "R2: 0.16\nC-index: 0.81\nLR chi2: 121",
                             "R2: 0.18\nC-index: 0.83\nLR chi2: 135",
                             "Percent of new information\nfrom including Adverse events: 11%")  ,
                  size = c(4, 4, 3.5, 3.5, 3.75), hjust = 0)

#The calibration plot displayes the bias corrected predictions versus the actual rates. 
#For patients with low risk of recurrance the model predict a low risk, but for patients 
#with higher risk, the predicted risk is overemphasized.             
#By adding adverse events to the model the model the explanatory power of the model 
#is increased by 11%. 

fit <- fit_lasso(X_, Y_, nfolds = 10, rule = "lambda.1se", seed = 123)
X2_ <- glmnet::makeX(data.frame(    
  TRTResection_w.o_wash=0, 
  TRTResection_w_wash= 1, 
  TRTHartmann_w.o_wash=0, 
  TRTHartmann_w_wash= 0, 
  TRTAPE.ELAPE=0,  
  Adverse.Events_Yes=0,  
  Radical.surgery_No.vs.Yes=0, 
  Neoadj.treatment_Yes.vs.No=1,
  Nodes_Yes.vs.No =0,
  Tumor.stage_3to4.vs.0to2 =1 ))
X3_ <- glmnet::makeX(data.frame(    
  TRTResection_w.o_wash=0, 
  TRTResection_w_wash= 1, 
  TRTHartmann_w.o_wash=0, 
  TRTHartmann_w_wash= 0, 
  TRTAPE.ELAPE=0,  
  Adverse.Events_Yes=1,  
  Radical.surgery_No.vs.Yes=0, 
  Neoadj.treatment_Yes.vs.No=1,
  Nodes_Yes.vs.No =0,
  Tumor.stage_3to4.vs.0to2 =1 ))
X4_ <- glmnet::makeX(data.frame(    
  TRTResection_w.o_wash=0, 
  TRTResection_w_wash= 1, 
  TRTHartmann_w.o_wash=0, 
  TRTHartmann_w_wash= 0, 
  TRTAPE.ELAPE=0,  
  Adverse.Events_Yes=0,  
  Radical.surgery_No.vs.Yes=1, 
  Neoadj.treatment_Yes.vs.No=1,
  Nodes_Yes.vs.No =0,
  Tumor.stage_3to4.vs.0to2 =1 ))
p1<-1- predict(fit, X_, Y_,  newx=X2_, pred.at=5)
p2<-1- predict(fit, X_, Y_,  newx=X3_, pred.at=5)
p3<-1- predict(fit, X_, Y_,  newx=X4_, pred.at=5)      

# For a patient without an adverse event the relative risk to an average patient is 
p1
p2
p3
#______________________________________________________________________________________________