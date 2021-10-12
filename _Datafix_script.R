rm(list = ls(all.names = TRUE))
gc()

library(tidyverse)
library(readxl)
library(lubridate)


path <- c("<insert path here>")    # path
A <- readRDS(paste0(path, "<Name of data set>.rds"))    # Name of data set


 A2 <- A %>%
   mutate(LR_dat = as_date(RDAT) , Journal_dat = as_date(Journalgrankning) ,
          Op_dat = as_date(Op.datum )      , Avliden_dat = as_date(DODSDATUM)) %>%
   mutate(tREC = replace_na(as.duration(Op_dat %--% LR_dat) / dyears(1) , 666)  ,
   tJOUR = replace_na(as.duration(Op_dat %--% Journal_dat) / dyears(1), 666 ) ,
   tDEATH = replace_na(as.duration(Op_dat %--% Avliden_dat) / dyears(1), 666) ) %>%
                   mutate( tREC = if_else(Patient.ID == "RAE1855", 666, tREC )) %>%
      rowwise() %>%
    mutate(FU = min(tREC, tJOUR, tDEATH)) %>%
    mutate(EVENT = case_when(  FU==tREC ~ 0 ,
                               FU==tJOUR ~ 1 ,
                               FU==tDEATH ~ 2) ) %>%
    mutate(Cens = case_when(EVENT == 0 ~ 0, EVENT != 0 ~ 1)) %>%
    mutate(Op.typ=na_if( Op.typ, "2")  ,
           pT = case_when(pT %in% c(0, 1, 2)  ~ 1 , pT %in% c(3, 4)  ~ 2) ,
           N.PAD = case_when(N.PAD == 0  ~ 1 , N.PAD %in% c(1, 2)  ~ 2) ,
           MKR=case_when(Mirko.Rad %in% c(0,2)~1, Mirko.Rad==1~2)) %>%
    dplyr::select(-Mirko.Rad) %>%
    mutate( Perop.rektskölj  =  na_if( Perop.rektskölj , "999") ,
            pT =  na_if( pT, "999") ,
            N.PAD =  na_if( N.PAD, "999") )    %>%

    mutate(SKOLJ = case_when(  Op.typ==1 &  Perop.rektskölj==0 ~ 0 ,
                               Op.typ==1 &  Perop.rektskölj==1 ~ 1 ,
                               Op.typ==4 &  Perop.rektskölj==0 ~ 0 ,
                               Op.typ==4 &  Perop.rektskölj==1 ~ 1 ,
                               Op.typ==3 ~ 2 ,
                               Op.typ==5 ~ 2    ) )   %>%
             mutate( Perop.rektskölj = factor(Perop.rektskölj) ,
               Adverse.Events = factor(Adverse.Events  )  ,    # , labels = c("Nej", "Ja")
            MKR = factor(MKR , ),                 #   labels = c("Radikalitet Nej", "Radikalitet Ja")
            N.PAD = factor(N.PAD) ,
            Neoadjuv.beh = factor(Neoadjuv.beh) ,
            Op.typ_ = factor( Op.typ  ) ,            #   , labels =  c("Resektion" , "APE/ELAPE", "Hartmann" ,  "Övrigt" )
            Recurrence = factor(Cens  )  ,    # , labels = c("Recurrence: No" , "Recurrence:  Yes")
            SKOLJ = factor(SKOLJ   ) ,                    #    labels = c("No" , "Yes" , "Not applicable")
            pT = factor(pT) )


  A3 <- A2 %>% dplyr::select( -c(Journalgrankning , Op.datum , LR.datum , Avliden.datum, DODSDATUM,
                            tREC, tJOUR, tDEATH ) )



