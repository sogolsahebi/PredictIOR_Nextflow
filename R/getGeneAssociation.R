##########################################################################
##########################################################################
## Get gene association (as continuous) with survival outcome (OS/PFS)
##########################################################################
##########################################################################

geneSurvCont <- function(dat.icb, clin = NULL, time.censor, missing.perc, const.int=0.001,
                         n.cutoff, feature, study, surv.outcome, cancer.type, treatment){

  if( !class(dat.icb) %in% c("SummarizedExperiment", "MultiAssayExperiment", "data.frame", "matrix") ){
    stop(message("function requires SummarizedExperiment, MultiAssayExperiment, data.frame, or matrix class of data"))
  }

  if( class(dat.icb) == "MultiAssayExperiment"){
    dat <- createSE(dat.icb)
    dat_expr <- assay(dat)
    dat_clin <- colData(dat)
  }

  if( class(dat.icb) == "SummarizedExperiment"){
    dat_expr <- assay(dat.icb)
    dat_clin <- colData(dat.icb)
  }

  if( sum(nrow(clin)) > 0  ){

    dat_expr <- dat.icb
    dat_clin <- clin

  }

        #cancer_type <- names( table( dat_clin$cancer_type )[ table( dat_clin$cancer_type ) >= n.cutoff ] )

        #message(paste(study))

        #data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type ]
        data <- dat_expr
        remove <- rem(data, missing.perc, const.int)

        if( length(remove) ){
          data <- data[-remove,]
        }

        data <- as.matrix( data[ rownames(data) %in% feature , ] )

        ## association with OS
        if( surv.outcome == "OS" ){

          if( nrow(data) & length( dat_clin$event_occurred_os[ !is.na( dat_clin$event_occurred_os ) ] ) >= n.cutoff ){

            res <- lapply(1:nrow(data), function(k){

              g <- as.numeric( scale( data[k , ] ) )
              names( g ) <- colnames( data )

              cox <- survCont( surv = dat_clin$event_occurred_os ,
                               time = dat_clin$survival_time_os ,
                               time.censor = time.censor ,
                               var = g )

              data.frame( Outcome = "OS",
                          Gene = rownames(data)[k],
                          Study = study,
                          Coef = cox["HR"],
                          SE = cox["SE"],
                          N = cox["N"],
                          Pval = cox["Pval"],
                          Cancer_type = cancer.type,
                          Treatment = treatment)

            })

            res <- do.call(rbind, res)
            rownames(res) <- NULL
            # res$FDR <- p.adjust(res$Pval, method = "BH")

          }else{  res <- data.frame( Outcome = "OS",
                                     Gene = NA,
                                     Study = study,
                                     Coef = NA,
                                     SE = NA,
                                     N = NA,
                                     Pval = NA,
                                     Cancer_type = NA,
                                     Treatment = NA)

          message("lack of number of samples and/or genes with known immunotherapy survival outcome")

          }

        }

        ## association with PFS
        if( surv.outcome == "PFS"){


          if( nrow(data) & length( dat_clin$event_occurred_pfs[ !is.na( dat_clin$event_occurred_pfs ) ] ) >= n.cutoff ){

            res <- lapply(1:nrow(data), function(k){

              g <- as.numeric( scale( data[k , ] ) )
              names( g ) <- colnames( data )

              cox <- survCont( surv = dat_clin$event_occurred_pfs ,
                               time = dat_clin$survival_time_pfs ,
                               time.censor = time.censor ,
                               var = g )

              data.frame( Outcome = "PFS",
                          Gene = rownames(data)[k],
                          Study = study,
                          Coef = cox["HR"],
                          SE = cox["SE"],
                          N = cox["N"],
                          Pval = cox["Pval"],
                          Cancer_type = cancer.type,
                          Treatment = treatment)

            })

            res <- do.call(rbind, res)
            rownames(res) <- NULL
           # res$FDR <- p.adjust(res$Pval, method = "BH")

          }else{  res <- data.frame( Outcome = "PFS",
                                     Gene = NA,
                                     Study = study,
                                     Coef = NA,
                                     SE = NA,
                                     N = NA,
                                     Pval = NA,
                                     Cancer_type = NA,
                                     Treatment = NA)

          message("lack of number of samples and/or genes with known immunotherapy survival outcome")

          }

        }

  return(res)
}


##########################################################################
##########################################################################
## Get gene association (as continuous) with survival outcome (OS/PFS)
##########################################################################
##########################################################################

geneSurvDicho <- function(dat.icb, clin = NULL, time.censor, missing.perc, const.int=0.001, n.cutoff, feature,
                          study, surv.outcome, n0.cutoff, n1.cutoff, method = "median", var.type, cancer.type, treatment){

  if( !class(dat.icb) %in% c("SummarizedExperiment", "MultiAssayExperiment", "matrix", "data.frame") ){

    stop(message("function requires SummarizedExperiment, MultiAssayExperiment, data.frame, or matrix class of data"))

  }

  if( class(dat.icb) == "MultiAssayExperiment"){


    dat <- createSE(dat.icb)
    dat_expr <- assay(dat)
    dat_clin <- colData(dat)

  }

  if( class(dat.icb) == "SummarizedExperiment"){

    dat_expr <- assay(dat.icb)
    dat_clin <- colData(dat.icb)

  }

  if( sum(nrow(clin)) > 0  ){

    dat_expr <- dat.icb
    dat_clin <- clin

  }

  #cancer_type <- names( table( dat_clin$cancer_type )[ table( dat_clin$cancer_type ) >= n.cutoff ] )

  #message(paste(study, cancer_type, sep="/"))

  #data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type]
  data <- dat_expr
  remove <- rem(data, missing.perc, const.int)

  if( length(remove) ){
    data <- data[-remove,]
  }

  data <- as.matrix( data[ rownames(data) %in% feature , ] )

  ## association with OS
  if( surv.outcome == "OS"){

    if( nrow(data) &  length( dat_clin$event_occurred_os[ !is.na( dat_clin$event_occurred_os ) ] ) >= n.cutoff ){

      res <- lapply(1:nrow(data), function(k){

        g <- as.numeric( scale( data[k , ] ) )
        names( g ) <- colnames( data )

        cox <- survDicho( surv = dat_clin$event_occurred_os ,
                          time = dat_clin$survival_time_os ,
                          time.censor = time.censor ,
                          var = g,
                          n0.cutoff = n0.cutoff,
                          n1.cutoff = n1.cutoff,
                          method = method,
                          var.type = var.type)

        data.frame( Outcome = "OS",
                    Gene = rownames(data)[k],
                    Study = study,
                    Coef = cox["HR"],
                    SE = cox["SE"],
                    N = cox["N"],
                    Pval = cox["Pval"],
                    Cancer_type = cancer.type,
                    Treatment = treatment)

      })

       res <- do.call(rbind, res)
       rownames(res) <- NULL
     #  res$FDR <- p.adjust(res$Pval, method = "BH")

     }else{  res <- data.frame( Outcome = "OS",
                               Gene = NA,
                               Study = study,
                               Coef = NA,
                               SE = NA,
                               N = NA,
                               Pval = NA,
                               Cancer_type = NA,
                               Treatment = NA)

    message("lack of number of samples and/or genes with known immunotherapy survival outcome")

    }

  }

  ## association with PFS
  if( surv.outcome == "PFS"){


    if( nrow(data) & length( dat_clin$event_occurred_pfs[ !is.na( dat_clin$event_occurred_pfs  ) ] ) >= n.cutoff ){

      res <- lapply(1:nrow(data), function(k){

        g <- as.numeric( scale( data[k , ] ) )
        names( g ) <- colnames( data )

        cox <- survDicho( surv = dat_clin$event_occurred_pfs ,
                          time = dat_clin$survival_time_pfs ,
                          time.censor= time.censor ,
                          var = g,
                          n0.cutoff = n0.cutoff,
                          n1.cutoff = n1.cutoff,
                          method = method,
                          var.type = var.type)

        data.frame( Outcome = "PFS",
                    Gene = rownames(data)[k],
                    Study = study,
                    Coef = cox["HR"],
                    SE = cox["SE"],
                    N = cox["N"],
                    Pval = cox["Pval"],
                    Cancer_type = cancer.type,
                    Treatment = treatment)

      })

      res <- do.call(rbind, res)
      rownames(res) <- NULL
     # res$FDR <- p.adjust(res$Pval, method = "BH")

    }else{  res <- data.frame( Outcome = "PFS",
                               Gene = NA,
                               Study = study,
                               Coef = NA,
                               SE = NA,
                               N = NA,
                               Pval = NA,
                               Cancer_type = NA,
                               Treatment = NA)

    message("lack of number of samples and/or genes with known immunotherapy survival outcome")

    }

  }

  return(res)
}

#################################################################
#################################################################
## Get gene association (as continuous) with response (R vs NR)
#################################################################
#################################################################
# n1.cutoff: cutoff for NR (== 1) samples
# n0.cutoff: cutoff for R (== 0) samples

geneLogReg <- function(dat.icb, clin = NULL, missing.perc, const.int=0.001, n.cutoff, feature, study,
                       n0.cutoff, n1.cutoff, cancer.type, treatment){

  if( !class(dat.icb) %in% c("SummarizedExperiment", "MultiAssayExperiment", "matrix", "data.frame") ){

    stop(message("function requires SummarizedExperiment, MultiAssayExperiment, data.frame, or matrix class of data"))

  }

  if( class(dat.icb) == "MultiAssayExperiment"){

    dat <- createSE(dat.icb)
    dat_expr <- assay(dat)
    dat_clin <- colData(dat)

  }

  if( class(dat.icb) == "SummarizedExperiment"){

    dat_expr <- assay(dat.icb)
    dat_clin <- colData(dat.icb)

  }

  if( sum(nrow(clin)) > 0  ){

    dat_expr <- dat.icb
    dat_clin <- clin

  }

    #cancer_type <- names( table( dat_clin$cancer_type )[ table( dat_clin$cancer_type ) >= n.cutoff ] )

    #message(paste(study, cancer_type, sep="/"))

    #data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type]
    data <- dat_expr
    remove <- rem(data, missing.perc, const.int)

    if( length(remove) ){
      data <- data[-remove,]
    }

    data <- as.matrix( data[ rownames(data) %in% feature , ] )

    if( nrow(data) & length(dat_clin$response) >= n.cutoff &
        sum(dat_clin$response == "NR", na.rm = TRUE) >= n1.cutoff &
        sum(dat_clin$response == "R", na.rm = TRUE) >= n0.cutoff ){

      res <- lapply(1:nrow(data), function(k){

        g <- as.numeric( scale( data[k , ] ) )
        names( g ) <- colnames( data )

        x <- ifelse( dat_clin$response %in% "R" , 0 ,
                     ifelse( dat_clin$response %in% "NR" , 1 , NA ) )

          fit <- glm( x ~ g , family=binomial( link="logit" ) )

          data.frame( Outcome = "R vs NR",
                      Gene = rownames(data)[k],
                      Study = study,
                      Coef = round( summary(fit)$coefficients[ "g" , "Estimate"  ] , 3 ),
                      SE = round( summary(fit)$coefficients[ "g" , "Std. Error" ] , 3 ),
                      N = length(x[!is.na(x)]),
                      Pval = summary(fit)$coefficients[ "g" , "Pr(>|z|)" ],
                      Cancer_type = cancer.type,
                      Treatment = treatment)

      })

      res <- do.call(rbind, res)
      rownames(res) <- NULL
    #  res$FDR <- p.adjust(res$Pval, method = "BH")

    }else{  res <- data.frame( Outcome = "R vs NR",
                               Gene = NA,
                               Study = study,
                               Coef = NA,
                               SE = NA,
                               N = NA,
                               Pval = NA,
                               Cancer_type = NA,
                               Treatment = NA)

    message("lack of number of samples and/or genes with known immunotherapy survival outcome")

    }

  return(res)

}

