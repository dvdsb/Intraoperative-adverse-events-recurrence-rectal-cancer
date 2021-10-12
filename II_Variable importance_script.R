
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

path <- c("<insert path here>")    # path
A <- readRDS(paste0(path, "<Name of data set>.rds"))    # Name of data set
B <- A %>% mutate(Recurrence = case_when(  EVENT2==1 ~ 1 , EVENT2==0 ~ 0 , EVENT2==2 ~ 0 ) ,
                  Recurrence = factor(Recurrence)) %>%
  dplyr::select(FU, Recurrence , Adverse.Events , MKR , N.PAD , Op.typ ,
                Neoadjuv.beh , Op.typ , SKOLJ , pT , Perop.rektskölj )


B <- B %>% mutate(Op.typ2 = case_when( Op.typ=="Resektion" &  Perop.rektskölj==0 ~ 1 ,   # Resektion utan skölj
                                       Op.typ=="Resektion" & Perop.rektskölj==1 ~ 2 ,   # Resektion med skölj
                                       Op.typ=="Hartmann" &  Perop.rektskölj==0 ~ 3 ,   # Hartmann utan skölj
                                       Op.typ=="Hartmann" & Perop.rektskölj==1 ~ 4 ,   # Hartmann med skölj
                                       Op.typ=="APE/ELAPE"  ~ 5 ))    %>%
  mutate( Op.typ_ = factor( Op.typ2 ))


REC <- recipe( Recurrence ~ . , data = B) %>%
  step_impute_knn( Adverse.Events , MKR , N.PAD , Neoadjuv.beh , Op.typ , Op.typ_ , pT ) %>%
  step_dummy( Adverse.Events , MKR , N.PAD ,  Neoadjuv.beh  , pT )

prepped_recipe <- prep(REC )
Bprep  <- juice(prepped_recipe)
names(Bprep)


C <- REC %>%
  # apply the recipe to the training data
  prep(B) %>%
  # extract the pre-processed training dataset
  juice()
C <- C %>%
mutate(  APE.ELAPE.vs.Resection = case_when( Op.typ=="APE/ELAPE" ~ 1,
                                             Op.typ=="Resektion"  ~ 0) ,
         APE.ELAPE.vs.Hartmann = case_when( Op.typ=="APE/ELAPE" ~ 1,
                                            Op.typ=="Hartmann"  ~ 0) ,
         Resection.vs.Hartmann = case_when( Op.typ=="Resektion" ~ 1,
                                            Op.typ=="Hartmann"  ~ 0) ,
         Resection_Wash.vs.No_Wash = case_when(  Op.typ_ == 2 ~ 1,
                                                 Op.typ_ == 1  ~ 0) ,
         Hartmann_Wash.vs.No_Wash = case_when(  Op.typ_ == 4 ~ 1,
                                                Op.typ_ == 3  ~ 0) )

#____________________________________________________________________________________________________________

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
                    APE.ELAPE.vs.Resection,
                    APE.ELAPE.vs.Hartmann ,
                    Resection.vs.Hartmann,
                    Resection_Wash.vs.No_Wash,
                    Hartmann_Wash.vs.No_Wash ,
                    Nodes_Yes.vs.No ,
                    Tumor.stage_3to4.vs.0to2 )) %>%
  mutate(EVENT = case_when( Recurrence == 0 ~ 0 , Recurrence == 1  ~ 1   )) %>%
  mutate(EVENT = factor(EVENT)) %>% dplyr::select(-c(Recurrence))

C2 <- tibble(time=survival::Surv(C$FU, C$EVENT, typ = "right")[,1],
             status=survival::Surv(C$FU, C$EVENT, typ = "right")[,2]) %>%
  cbind(C)

# A. APE.ELAPE.vs.Resection,
# B. APE.ELAPE.vs.Hartmann ,
# C. Resection.vs.Hartmann,
# D. Resection_Wash.vs.No_Wash,
# E. Hartmann_Wash.vs.No_Wash

C2_0 <- C2 %>% mutate(TRT = Op.typ_) %>%
      dplyr::select(-c(APE.ELAPE.vs.Resection,
                       APE.ELAPE.vs.Hartmann ,
                       Resection.vs.Hartmann,
                       Resection_Wash.vs.No_Wash,
                       Hartmann_Wash.vs.No_Wash , Op.typ_ ))

C2_A <- C2 %>% filter(APE.ELAPE.vs.Resection %in% c(0, 1)) %>%
    mutate(TRT = APE.ELAPE.vs.Resection) %>%       dplyr::select(-c(APE.ELAPE.vs.Resection,
                                                                    APE.ELAPE.vs.Hartmann ,
                                                                    Resection.vs.Hartmann,
                                                                    Resection_Wash.vs.No_Wash,
                                                                    Hartmann_Wash.vs.No_Wash , Op.typ_ ))
C2_B <- C2 %>% filter(APE.ELAPE.vs.Hartmann %in% c(0, 1)) %>%
    mutate(TRT = APE.ELAPE.vs.Hartmann) %>%       dplyr::select(-c(APE.ELAPE.vs.Resection,
                                                                   APE.ELAPE.vs.Hartmann ,
                                                                   Resection.vs.Hartmann,
                                                                   Resection_Wash.vs.No_Wash,
                                                                   Hartmann_Wash.vs.No_Wash , Op.typ_))
C2_C <- C2 %>% filter(Resection.vs.Hartmann %in% c(0, 1)) %>%
    mutate(TRT = Resection.vs.Hartmann) %>%       dplyr::select(-c(APE.ELAPE.vs.Resection,
                                                                   APE.ELAPE.vs.Hartmann ,
                                                                   Resection.vs.Hartmann,
                                                                   Resection_Wash.vs.No_Wash,
                                                                   Hartmann_Wash.vs.No_Wash , Op.typ_ ))
C2_D <- C2 %>% filter(Resection_Wash.vs.No_Wash %in% c(0, 1)) %>%
    mutate(TRT = Resection_Wash.vs.No_Wash) %>%       dplyr::select(-c(APE.ELAPE.vs.Resection,
                                                                       APE.ELAPE.vs.Hartmann ,
                                                                       Resection.vs.Hartmann,
                                                                       Resection_Wash.vs.No_Wash,
                                                                       Hartmann_Wash.vs.No_Wash, Op.typ_ ))
C2_E <- C2 %>% filter(Hartmann_Wash.vs.No_Wash %in% c(0, 1)) %>%
    mutate(TRT = Hartmann_Wash.vs.No_Wash) %>%       dplyr::select(-c(APE.ELAPE.vs.Resection,
                                                                      APE.ELAPE.vs.Hartmann ,
                                                                      Resection.vs.Hartmann,
                                                                      Resection_Wash.vs.No_Wash,
                                                                      Hartmann_Wash.vs.No_Wash , Op.typ_ ))

# 1. Determine lambda.1se by cross-validation
# 2. Store which variables selected when estimating with lambda.1se
# 3. Take the selected variables and estimate a standard Cox regression.

HD_res <- function(splits ) {

  mod <- cv.glmnet(   survival::Surv(time, status, typ = "right") ~
                        Adverse.Events_Yes +
                        Radical.surgery_No.vs.Yes +
                        Neoadj.treatment_Yes.vs.No +
                        TRT +
                        Nodes_Yes.vs.No +
                        Tumor.stage_3to4.vs.0to2  ,
                      family="cox" , alpha = 1  , nfolds = 10 , se.model.frame = TRUE ,
                      type.measure="C" ,
                      data = analysis(splits)  )

  crs <- coef(mod , s = mod$lambda.1se)                      #    mod$lambda.1se
  VRS <- dimnames(crs)[[1]][which(crs != 0)]

  hd_ <- analysis(splits)

  CXRegDat <- hd_ %>% dplyr::select(  any_of(c("time", "status", VRS)) )

  GLM_ <-  coxph(survival::Surv(time, status, typ = "right") ~ .,
            data=CXRegDat , x = TRUE, y=TRUE )
}

#__________________________________________________________________________________________________________
# 1. Do the previous steps for each bootstrap replicate
# 2. Store which variables were selected and their standard HR´s (except for the intercept)
HD_res2 <-  function( dat_ , Ntimes )
  {
  set.seed(3352)
    N_  <- Ntimes
  btsrp <- bootstraps( dat_  , times = N_ , strata = EVENT )

  OUT <- btsrp %>% mutate(MDL = map(splits, HD_res ) ) %>%
    mutate(coef_info = map(MDL, tidy) )

  OUT2 <- OUT %>%
    unnest(coef_info)

  OUT3 <- tibble(OUT2$id , OUT2$term , OUT2$estimate , OUT2$p.value)

  OUT4 <- OUT3  %>% dplyr::filter ( `OUT2$term`  != "(Intercept)")
  return(OUT4)
  }

SUM_C2_0x <- HD_res2(dat_ = C2_0 , Ntimes = 1000 )
SUM_C2_Ax <- HD_res2(dat_ = C2_A , Ntimes = 1000 )
SUM_C2_Bx <- HD_res2(dat_ = C2_B , Ntimes = 1000 )
SUM_C2_Cx <- HD_res2(dat_ = C2_C , Ntimes = 1000 )
SUM_C2_Dx <- HD_res2(dat_ = C2_D , Ntimes = 1000 )
SUM_C2_Ex <- HD_res2(dat_ = C2_E , Ntimes = 1000 )

                                                    #__________________________________
                                          # A. APE.ELAPE.vs.Resection,
                                          # B. APE.ELAPE.vs.Hartmann ,
                                          # C. Resection.vs.Hartmann,
                                          # D. Resection_Wash.vs.No_Wash,
                                          # E. Hartmann_Wash.vs.No_Wash
                                                        table(SUM_C2_0x$`OUT2$term`)
                                                        table(SUM_C2_Ax$`OUT2$term`)
                                                        table(SUM_C2_Bx$`OUT2$term`)
                                                        table(SUM_C2_Cx$`OUT2$term`)
                                                        table(SUM_C2_Dx$`OUT2$term`)
                                                        table(SUM_C2_Ex$`OUT2$term`)
                                                    #__________________________________

# Make the data in another format
HD_res22 <-  function( dat_ , a )
{     DAT <- dat_
  OUT_ <- tibble(Var1 = DAT[2] ,  EST =  DAT[3], A = a)  %>%
         filter(Var1 != "TRT")
 return(OUT_)  }

    HD_res222 <-  function( dat_ , a )
    {     DAT <- dat_
    OUT_ <- tibble(Var1 = DAT[2] ,  EST =  DAT[3], A = a)  %>%
      filter(Var1 == "TRT")
    return(OUT_)  }

SUM_C2_0xx <- HD_res22(dat_ = SUM_C2_0x , a = 0 )
SUM_C2_Axx <- HD_res222(dat_ = SUM_C2_Ax , a = 1 )
SUM_C2_Bxx <- HD_res222(dat_ = SUM_C2_Bx , a = 2 )
SUM_C2_Cxx <- HD_res222(dat_ = SUM_C2_Cx , a = 3 )
SUM_C2_Dxx <- HD_res222(dat_ = SUM_C2_Dx , a = 4 )
SUM_C2_Exx <- HD_res222(dat_ = SUM_C2_Ex , a = 5 )

# Combine it and save.
SUM_C3 <- rbind(SUM_C2_0xx ,
             SUM_C2_Axx ,
             SUM_C2_Bxx ,
             SUM_C2_Cxx ,
             SUM_C2_Dxx ,
             SUM_C2_Exx)

SUM_C3 <- SUM_C3 %>%
  mutate(V = case_when(Var1 =="Adverse.Events_Yes" ~ "Adverse.Events_Yes" ,
                       Var1 =="Neoadj.treatment_Yes.vs.No" ~ "Neoadj.treatment_Yes.vs.No" ,
                       Var1 =="Nodes_Yes.vs.No" ~ "Nodes_Yes.vs.No" ,
                       Var1 =="Radical.surgery_No.vs.Yes" ~ "Radical.surgery_No.vs.Yes" ,
                       Var1 =="Tumor.stage_3to4.vs.0to2" ~ "Tumor.stage_3to4.vs.0to2 " ,
                       A ==1 ~ "APE.ELAPE.vs.Resection" ,
                       A ==2 ~ "APE.ELAPE.vs.Hartmann" ,
                       A ==3 ~ "Resection.vs.Hartmann" ,
                       A ==4 ~ "Resection_Wash.vs.No_Wash" ,
                       A ==5 ~ "Hartmann_Wash.vs.No_Wash" ))


#______________________________________________________________________________________________________

SUM_C3 <- SUM_C3 %>%
            dplyr::select(-Var1 )
SUM_C3 <- SUM_C3 %>%
  mutate(ordr = case_when( V == "Adverse.Events_Yes"    ~  1,
                           V == "Radical.surgery_No.vs.Yes"  ~  2,
                           V == "Tumor.stage_3to4.vs.0to2 "  ~  3,
                           V == "Nodes_Yes.vs.No"  ~  4,
                           V == "Neoadj.treatment_Yes.vs.No"  ~  5,
                           V == "APE.ELAPE.vs.Resection"  ~  6,
                           V == "APE.ELAPE.vs.Hartmann"  ~  7,
                           V == "Resection.vs.Hartmann"  ~  8,
                           V == "Resection_Wash.vs.No_Wash"  ~  9,
                           V == "Hartmann_Wash.vs.No_Wash"  ~ 10 ),
         V_ = factor(ordr ,  labels = c("Adverse events Yes vs No",
                                        "Radical surgery No vs Yes",
                                        "Tumour stage T3-4 vs T0-2",
                                        "Lymph node metastasis Yes vs No",
                                        "Neoadjuvant treatment Yes vs No",
                                        "APE vs AR",
                                        "APE vs HP",
                                        "AR vs HP",
                                        "AR Washout vs  No washout" ,
                                        "HP Washout vs No washout"  ) ),
         HR = exp(EST$`OUT2$estimate`),
         logHR = EST$`OUT2$estimate`)

SUMMARY <- SUM_C3 %>%
  group_by(V, V_ , ordr) %>%
  summarise(median = quantile(HR, na.rm=T, probs = 0.5),
            q1 = quantile(HR, na.rm=T, probs = 0.25),
            q3 = quantile(HR, na.rm=T, probs = 0.75),
            n=100*(n()/1000)) %>%
            mutate(res=paste0(signif(median, digits = 2),
                  "(",signif(q1, digits = 2),";",signif(q3, digits = 2),")"),
                  res2=paste0(signif(n, digits = 2), "%"))

SUM_C3$V_x <-  fct_reorder(SUM_C3$V_, SUM_C3$ordr, min)
SUMMARY$V_x <-  fct_reorder(SUMMARY$V_, SUMMARY$ordr, min)



#1. Adverse events Yes vs No
#2. Radical surgery No vs Yes
#3. Tumour stage T3-4 vs T0-2
#4. Lymph node metastasis Yes vs No
#5. Neoadjuvant treatment Yes vs No
#6. APE vs AR
#7. APE vs HP
#8. AR vs HP
#9: AR Washout vs  No washout
#10: HP Washout vs No washout

TExTT <- SUMMARY %>%
  mutate(Result_ = paste(res2, " \n ",res, sep='')) %>%
  dplyr::select(V, ordr, Result_  ) %>%
  arrange(ordr)

TExTT

 #_____________________________________________________________________________________________

tag_facet <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf,
                      hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust,
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE) }

dev.new(width=20, height=20)
window()
pd <- position_dodge(.25)
par(mar = c(0.4, 0.1, 0.4, 0.1))


PLOTT <- SUM_C3  %>% ggplot( aes(x=logHR)) +
  geom_histogram(alpha=0.5, bins=100) +
  theme_bw() +
  geom_rug() +
  facet_wrap( ~  V_x , nrow = 5 )  +
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text.x = element_text(size = 7, colour = "black", angle = 0)) +
  xlab("log(Hazard ratio)") + ylab("Frequency") +
   scale_y_continuous(limits=c(0 , 50)) +
  scale_x_continuous(limits=c(-2 , 2.5)) +
  geom_vline(xintercept = 0, linetype=2,
             color = "grey")

my_tag <- TExTT[4]    # 3
my_tag_t <- t(my_tag)

PLOTT2 <- tag_facet(PLOTT,
                    x = -1.5, y = 25,
                    vjust = -1, hjust = 0.25,
                    open = "", close = "",
                    fontface = 1,
                    size = 2.55,
                    family = "arial",
                    tag_pool = my_tag_t)


PLOTT2


####################################################################################################################
