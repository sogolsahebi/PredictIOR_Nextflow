## load libraries

library(meta)
library(metafor)
library(forestplot)
library(ggplot2)
library(ggrepel)

##############################################################
##############################################################
## create a function to specify cancer specific analyses
##############################################################
##############################################################

cancer.mod <- function( cancer ){

  cancer$Coef <- as.numeric( as.character( cancer$Coef ) )
  cancer$Study <- as.character( cancer$Study )
  cancer$Cancer_type <- as.character( cancer$Cancer_type )
  cancer$Pval <- as.numeric(as.character( cancer$Pval ))
  cancer$SE <- as.numeric(as.character( cancer$SE ))

  tab <- table( cancer$Cancer_type)[ table(cancer$Cancer_type) %in% c(1,2) ]
  cancer$Cancer_type[ cancer$Cancer_type %in% names(tab) ] <- "Other"
  cancer$Study <- as.character(paste(cancer$Study, ", n = ", cancer$N, sep=""))

  cancer

}

##############################################################
##############################################################
## create a function to specify treatment analyses
##############################################################
##############################################################

treatment.mod <- function( treatment ){

  treatment$Coef <- as.numeric( as.character( treatment$Coef ) )
  treatment$Study <- as.character( treatment$Study )
  treatment$Cancer_type <- as.character( treatment$Cancer_type )
  treatment$Pval <- as.numeric(as.character( treatment$Pval ))
  treatment$SE <- as.numeric(as.character( treatment$SE ))
  treatment$Study <- as.character(paste(treatment$Study, ", n = ", treatment$N, sep=""))

  tab <- table( treatment$Treatment)[ table(treatment$Treatment) %in% c(1,2) ]
  treatment$Treatment[ treatment$Treatment %in% names(tab) ] <- "Other"

  treatment

}

##############################################################
##############################################################
## create meta-analysis function
##############################################################
##############################################################

meta.fun <- function(coef, se, study, pval, n, cancer.type, treatment, cancer.spec = FALSE, treatment.spec = FALSE, feature){

  data <- data.frame( Gene = feature,
                      Study = as.character( study ),
                      N = n,
                      Coef = as.numeric(as.character( coef )),
                      SE = as.numeric(as.character( se )),
                      Pval = as.numeric(as.character( pval )),
                      Cancer_type = as.character( cancer.type ),
                      Treatment = as.character(treatment))

  data <- data[ order( data$Coef ) , ]

  ## at least 3 studies needed to do the random effect meta-analyses
  if(nrow(data) >= 3){

    if( cancer.spec == TRUE ){

      cancer <- cancer.mod( data )
      data <- cancer[ order( cancer$Coef ) , ]

    }

    if( treatment.spec == TRUE ){

      treatment <- treatment.mod( data )
      data <- treatment[ order( treatment$Coef ) , ]

    }

    meta <- metagen( TE = Coef,
                     seTE = SE ,
                     data = data ,
                     studlab = Study ,
                     fixed = FALSE ,
                     random = TRUE ,
                     control = list( maxiter = 10000 , stepadj=0.5 ) )

    meta_res <- data.frame(Gene = feature,
                           Coef = meta$TE.random ,
                           SE = meta$seTE.random ,
                           CI_lower = meta$lower.random ,
                           CI_upper = meta$upper.random ,
                           Pval = meta$pval.random ,
                           I2 = meta$I2 ,
                           Q_Pval = meta$pval.Q )
  }else{

    print("not enough studies to do meta-analysis")
    meta <- NA
    meta_res <- data.frame(Gene = feature,
                           Coef = NA,
                           SE =  NA,
                           CI_lower = NA,
                           CI_upper = NA,
                           Pval = NA,
                           I2= NA,
                           Q_Pval = NA)
         }

  return(list(input_data = data,
              meta_output = meta,
              meta_summery = meta_res))
  }

##############################################################
##############################################################
## create meta-analysis per-cancer function
##############################################################
##############################################################

metaPerCan.fun <- function(coef, se, study, pval, n, cancer.type, treatment, cancer.spec = TRUE, feature){

  data <- data.frame( Gene = feature,
                      Study = as.character( study ),
                      N = n,
                      Coef = as.numeric(as.character( coef )),
                      SE = as.numeric(as.character( se )),
                      Pval = as.numeric(as.character( pval )),
                      Cancer_type = as.character( cancer.type ),
                      Treatment = as.character(treatment))

  data <- data[ order( data$Coef ) , ]

  ## at least 3 studies needed to do the random effect meta-analyses
  if(nrow(data) >= 3){

    cancer <- cancer.mod( data )
    data <- cancer[ order( cancer$Coef ) , ]

    group <- names(table(data$Cancer_type)[table(data$Cancer_type) >= 3])

    if( length(group) > 0){

      res <- lapply(1:length(group), function(i){

        sub_data <- data[data$Cancer_type == group[i], ]

        meta <- metagen( TE = Coef,
                         seTE = SE ,
                         data = sub_data ,
                         studlab = sub_data$Study ,
                         fixed = FALSE ,
                         random = TRUE ,
                         control = list( maxiter = 10000 , stepadj=0.5 ) )

        meta_res <- data.frame(Cancer_type = group[i],
                               Gene = feature,
                               Coef = meta$TE.random ,
                               SE = meta$seTE.random ,
                               CI_lower = meta$lower.random ,
                               CI_upper = meta$upper.random ,
                               Pval = meta$pval.random ,
                               I2 = meta$I2 ,
                               Q_Pval = meta$pval.Q )

        list(input_data = sub_data,
             meta_output = meta,
             meta_summery = meta_res)

      })

      names(res) <- group

    }else{

      print("not enough studies to do cancer-specific meta-analysis")
      res <- NA

    }

  }else{

    print("not enough studies to do cancer-specific meta-analysis")
    res <- NA

  }

  return(res)

}

##############################################################
##############################################################
## create meta-analysis per-treatment function
##############################################################
##############################################################

metaPerTreatment.fun <- function(coef, se, study, pval, n, cancer.type, treatment, treatment.spec = TRUE, feature){

  data <- data.frame( Gene = feature,
                      Study = as.character( study ),
                      N = n,
                      Coef = as.numeric(as.character( coef )),
                      SE = as.numeric(as.character( se )),
                      Pval = as.numeric(as.character( pval )),
                      Cancer_type = as.character( cancer.type ),
                      Treatment = as.character(treatment))

  data <- data[ order( data$Coef ) , ]

  ## at least 3 studies needed to do the random effect meta-analyses
  if(nrow(data) >= 3){

    treatment <- treatment.mod( data )
    data <- treatment[ order( treatment$Coef ) , ]

    group <- names(table(data$Treatment)[table(data$Treatment) >= 3])

    if( length(group) > 0){

      res <- lapply(1:length(group), function(i){

        sub_data <- data[data$Treatment == group[i], ]

        meta <- metagen( TE = Coef,
                         seTE = SE ,
                         data = sub_data ,
                         studlab = sub_data$Study ,
                         fixed = FALSE ,
                         random = TRUE ,
                         control = list( maxiter = 10000 , stepadj=0.5 ) )

        meta_res <- data.frame(Treatment = group[i],
                               Gene = feature,
                               Coef = meta$TE.random ,
                               SE = meta$seTE.random ,
                               CI_lower = meta$lower.random ,
                               CI_upper = meta$upper.random ,
                               Pval = meta$pval.random ,
                               I2 = meta$I2 ,
                               Q_Pval = meta$pval.Q )

        list(input_data = sub_data,
             meta_output = meta,
             meta_summery = meta_res)

      })

      names(res) <- group

    }else{

      print("not enough studies to do cancer-specific meta-analysis")
      res <- NA

    }

  }else{

    print("not enough studies to do cancer-specific meta-analysis")
    res <- NA

  }

  return(res)

}


##############################################################
##############################################################
## create forestplot function: Pan cancer analysis
##############################################################
##############################################################

forestPlot <- function(coef, se, study, pval, n , cancer.type, treatment, xlab , label, feature){

  study <- paste( study , ", n = " , n , sep= "" )
  res <- meta.fun(coef, se, study, pval, n, cancer.type, treatment, cancer.spec = FALSE, treatment.spec = FALSE, feature)
  data <- res$input_data
  meta <- res$meta_output

  m <- c( min( c( 0 , data$Coef ) , na.rm=TRUE) - .5 , ( max( c( 0 , abs(data$Coef) ) , na.rm=TRUE) ) + .5 )

  forest( meta ,
          leftcols = c("studlab", "effect.ci" , "Pval" ),
          leftlabs= c( "Study" , paste(label, "[95%CI]", sep = " ") , "P-value" ) ,
          xlab = xlab,
          digits.se = 2 ,
          colgap.forest=unit(10, "mm") ,
          plotwidth = unit( 30 , "mm") ,
          pooled.totals = TRUE,
          smlab = " ",
          comb.random =TRUE,
          comb.fixed = FALSE,
          text.fixed.w = FALSE,
          layout = "JAMA",
          print.I2.ci = TRUE,
          print.Q = FALSE,
          print.pval.Q = TRUE,
          print.I2 = TRUE,
          print.tau2 = FALSE,
          resid.hetstat = FALSE,
          test.overall.random = TRUE,
          test.overall.fixed = FALSE,
          xlim = m ,
          col.square= "#4d4d4d" ,
          col.study= "#4d4d4d" ,
          col.square.lines = "#4d4d4d" ,
          col.diamond.random  = "#2166ac"  ,
          col.diamond.lines.random  ="#2166ac" ,
          col.by = "#2166ac",
          addrow.subgroups=TRUE )

  }

#########################################################
#########################################################
## create forestplot function: Per cancer analysis
#########################################################
#########################################################

forestPlotPerCan <- function( coef, se, study, pval, n, cancer.type, treatment, xlab , label, feature){

  res <- meta.fun(coef, se, study, pval, n, cancer.type, treatment, cancer.spec = TRUE, treatment.spec = FALSE, feature)
  cancer <- cancer.mod(res$input_data)

  remove <- names( table( cancer$Cancer_type)[ table(cancer$Cancer_type) %in% c(1,2) ] )

  if(length(table( cancer$Cancer_type )[ table(cancer$Cancer_type) >= 3 ]) > 1){

    if( length( unique( cancer$Cancer_type[ !cancer$Cancer_type %in% remove ] ) ) > 1 ){

      m <- c( min( c( 0 , cancer$Coef ) , na.rm=TRUE) - .5 , ( max( c( 0 , abs(cancer$Coef) ) , na.rm=TRUE) ) + .5 )
      meta <- res$meta_output

      if( length(remove) > 0 ){

        meta.subgroup <- update(meta ,
                                     byvar = Cancer_type ,
                                     exclude = cancer$cancer_type %in% remove ,
                                     fixed = FALSE ,
                                     random = TRUE ,
                                     control = list( maxiter = 10000 , stepadj=0.5 ) )
      } else{
        meta.subgroup <- update(meta ,
                                     byvar = Cancer_type ,
                                     comb.random = TRUE ,
                                     fixed = FALSE ,
                                     random = TRUE ,
                                     control = list( maxiter = 10000 , stepadj=0.5 ) )
      }

    }

  }else{

    stop("not enough studies to do cancer-specific sub-group meta-analysis")

  }

  forest( meta.subgroup,
          digits.se = 2,
          colgap.forest=unit(10, "mm") ,
          plotwidth = unit( 30 , "mm") ,
          xlab = xlab,
          pooled.totals = TRUE,
          smlab = " ",
          comb.random =TRUE,
          comb.fixed = FALSE,
          text.fixed.w = FALSE,
          layout = "JAMA",
          print.I2.ci = TRUE,
          print.Q = FALSE,
          print.pval.Q = TRUE,
          print.I2 = TRUE,
          print.tau2 = FALSE,
          resid.hetstat = FALSE,
          test.overall.random = TRUE,
          test.overall.fixed = FALSE,
          xlim = m ,
          col.square= "#4d4d4d" ,
          col.study= "#4d4d4d" ,
          col.square.lines = "#4d4d4d" ,
          col.diamond.random  = "#2166ac"  ,
          col.diamond.lines.random  ="#2166ac" ,
          col.by = "#2166ac",
          addrow.subgroups=TRUE )

}

#########################################################
#########################################################
## create forestplot function: Per treatment analysis
#########################################################
#########################################################

forestPlotPerTreatment <- function( coef, se, study, pval, n , cancer.type, treatment, feature, xlab , label){

  res <- meta.fun(coef, se, study, pval, n, cancer.type, treatment, cancer.spec = FALSE, treatment.spec = TRUE, feature)
  treatment <- treatment.mod(res$input_data)

  remove <- names( table( treatment$Treatment )[ table(treatment$Treatment) %in% c(1,2) ] )

  if(length(table( treatment$Treatment )[ table(treatment$Treatment) >= 3 ]) > 1){

    if( length( unique( treatment$Treatment[ !treatment$Treatment %in% remove ] ) ) > 1 ){

      m <- c( min( c( 0 , treatment$Coef ) , na.rm=TRUE) - .5 , ( max( c( 0 , abs(treatment$Coef) ) , na.rm=TRUE) ) + .5 )
      meta <- res$meta_output

      if( length(remove) > 0 ){

        meta.subgroup <- update(meta ,
                                     byvar = Treatment ,
                                     exclude = treatment$Treatment %in% remove ,
                                     fixed = FALSE ,
                                     random = TRUE ,
                                     control = list( maxiter = 10000 , stepadj=0.5 ) )
      } else{
        meta.subgroup <- update(meta ,
                                     byvar = Treatment ,
                                     comb.random = TRUE ,
                                     fixed = FALSE ,
                                     random = TRUE ,
                                     control = list( maxiter = 10000 , stepadj=0.5 ) )
      }

    }

  }else{

    stop("not enough studies to do treatment-specific sub-group meta-analysis")

  }

  forest( meta.subgroup,
          #leftcols = c("studlab", "effect.ci", "Pval"),
          #leftlabs= c( "Study" , paste(label, "[95%CI]", sep = " ") , "P-value" ) ,
          digits.se = 2,
          colgap.forest=unit(10, "mm") ,
          plotwidth = unit( 30 , "mm") ,
          xlab = xlab,
          pooled.totals = TRUE,
          smlab = " ",
          comb.random =TRUE,
          comb.fixed = FALSE,
          text.fixed.w = FALSE,
          layout = "JAMA",
          print.I2.ci = TRUE,
          print.Q = FALSE,
          print.pval.Q = TRUE,
          print.I2 = TRUE,
          print.tau2 = FALSE,
          resid.hetstat = FALSE,
          test.overall.random = TRUE,
          test.overall.fixed = FALSE,
          xlim = m ,
          col.square= "#4d4d4d" ,
          col.study= "#4d4d4d" ,
          col.square.lines = "#4d4d4d" ,
          col.diamond.random  = "#2166ac"  ,
          col.diamond.lines.random  ="#2166ac" ,
          col.by = "#2166ac",
          addrow.subgroups=TRUE )


 }



