## load libraries

library(GSVA)

##########################################################
##########################################################
# scale data
##########################################################
##########################################################

scale.fun <- function( x ){

  rid = rownames(x)
  cid = colnames(x)
  out = t( apply( x , 1 , scale ) )
  rownames(out) = rid
  colnames(out) = cid
  out

}

#####################################################################
## Get signature score: GSVA
#####################################################################

geneSigGSVA <- function(dat.icb, clin = NULL, sig, sig.name, missing.perc, const.int = 0.001, n.cutoff, sig.perc, study){

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

    #data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type]
    data <- dat_expr
    remove <- rem(data, missing.perc, const.int)

    if( length(remove) ){
      data <- data[-remove,]
    }

    if( ifelse( is.null( nrow( data[ rownames(data) %in% sig$gene_name , ]) ) , 1 , nrow( data[ rownames(data) %in% sig$gene_name , ] ) ) / length( sig$gene_name ) > sig.perc & ncol(data) >= n.cutoff ){

      #print( paste( signature_name , "|" , "GSVA" , sep=" " ) )

      geneSig <- NULL
      gsvaPar <- gsvaParam(scale.fun( x=data ) , list(sig$gene_name))
      geneSig <- gsva(gsvaPar, verbose=FALSE)
      #geneSig <- as.numeric(gsva( scale.fun( x=data ) , list(sig$gene_name) , verbose=FALSE ) )
      #names( geneSig ) <- colnames(data)

    }else{

      geneSig <- NA
      message("not enough samples and/or genes in a data")

    }

  return(geneSig)
}


#####################################################################
## Get signature score: ssGSEA
#####################################################################

geneSigssGSEA <- function(dat.icb, clin = NULL, sig, sig.name, missing.perc, const.int = 0.001, n.cutoff, sig.perc, study){

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
  #data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type]

  data <- dat_expr
  remove <- rem(data, missing.perc, const.int)

  if( length(remove) ){
    data <- data[-remove,]
  }

  if( ifelse( is.null( nrow( data[ rownames(data) %in% sig$gene_name , ]) ) , 1 , nrow( data[ rownames(data) %in% sig$gene_name , ] ) ) / length( sig$gene_name ) > sig.perc & ncol(data) >= n.cutoff ){

    #print( paste( signature_name , "|" , "ssGSEA" , sep=" " ) )

    geneSig <- NULL
    geneSig <- as.numeric(gsva( scale.fun( x=data ) , list(sig$gene_name) ,
                                method = "ssgsea", kcdf = "Gaussian", verbose=FALSE ) )
    names( geneSig ) <- colnames(data)

  }else{

    geneSig <- NA
    message("not enough samples and/or genes in a data")

  }

  return(geneSig)
}

#####################################################################
## Get signature score: (weighted) mean
#####################################################################

geneSigMean <- function(dat.icb, clin = NULL, sig, sig.name, missing.perc, const.int =0.001, n.cutoff, sig.perc, study){

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

    if( ifelse( is.null( nrow( data[ rownames(data) %in% sig$gene_name , ]) ) , 1 , nrow( data[ rownames(data) %in% sig$gene_name , ] ) ) / length( sig$gene_name ) > sig.perc & ncol(data) >= n.cutoff ){

      #print( paste( signature_name , "|" , "Weighted Mean" , sep=" " ) )

      gene <- intersect( rownames(data) , sig$gene_name)
      s <- sig[ sig$gene_name %in% gene, ]

      scaled_dat <- scale.fun( x= data[ gene , ] )

      remove <- which(is.na(scaled_dat))
      if(length(remove)){
        scaled_dat <- scaled_dat[-which(is.na(scaled_dat)), ]
      }else{ scaled_dat }

      geneSig <- NULL
      geneSig <- apply( scaled_dat , 2 , function(x) ( sum( x * s$weight, na.rm=TRUE  )  /  nrow( s ) ) )
      names( geneSig ) <- colnames(data)

    }else{

      geneSig <- NA
      message("not enough samples and/or genes in a data")

    }

  return(geneSig)
}


#####################################################################
## Get signature score: Median
#####################################################################

geneSigMedian <- function(dat.icb, clin = NULL, sig, sig.name, missing.perc, const.int =0.001, n.cutoff, sig.perc, study){

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

  if( ifelse( is.null( nrow( data[ rownames(data) %in% sig$gene_name , ]) ) , 1 , nrow( data[ rownames(data) %in% sig$gene_name , ] ) ) / length( sig$gene_name ) > sig.perc & ncol(data) >= n.cutoff ){

    #print( paste( signature_name , "|" , "Median" , sep=" " ) )

    gene <- intersect( rownames(data) , sig$gene_name)
    s <- sig[ sig$gene_name %in% gene, ]

    scaled_dat <- scale.fun( x= data[ gene , ] )

    remove <- which(is.na(scaled_dat))
    if(length(remove)){
      scaled_dat <- scaled_dat[-which(is.na(scaled_dat)), ]
    }else{ scaled_dat }

    geneSig <- NULL
    geneSig <- apply( scaled_dat , 2 , function(x) ( median(x, na.rm=TRUE) ) )
    names( geneSig ) <- colnames(data)

  }else{

    geneSig <- NA
    message("not enough samples and/or genes in a data")

  }

  return(geneSig)
}


#####################################################################
## Get signature score: (weighted) summation
#####################################################################

geneSigSum <- function(dat.icb, clin = NULL, sig, sig.name, missing.perc, const.int =0.001, n.cutoff, sig.perc, study){

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

  if( ifelse( is.null( nrow( data[ rownames(data) %in% sig$gene_name , ]) ) , 1 , nrow( data[ rownames(data) %in% sig$gene_name , ] ) ) / length( sig$gene_name ) > sig.perc & ncol(data) >= n.cutoff ){

    #print( paste( signature_name , "|" , "Median" , sep=" " ) )

    gene <- intersect( rownames(data) , sig$gene_name)
    s <- sig[ sig$gene_name %in% gene, ]

    scaled_dat <- scale.fun( x= data[ gene , ] )

    remove <- which(is.na(scaled_dat))
    if(length(remove)){
      scaled_dat <- scaled_dat[-which(is.na(scaled_dat)), ]
    }else{ scaled_dat }

    geneSig <- NULL
    geneSig <- apply( scaled_dat , 2 , function(x) ( sum( x * s$weight , na.rm=TRUE ) ) )
    names( geneSig ) <- colnames(data)

  }else{

    geneSig <- NA
    message("not enough samples and/or genes in a data")

  }

  return(geneSig)
}

#####################################################################
## Get signature score: Geometric mean
#####################################################################

geneSigGMean <- function(dat.icb, clin = NULL, sig, sig.name, missing.perc, const.int =0.001, n.cutoff, sig.perc, study){

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

  if( ifelse( is.null( nrow( data[ rownames(data) %in% sig$gene_name , ]) ) , 1 , nrow( data[ rownames(data) %in% sig$gene_name , ] ) ) / length( sig$gene_name ) > sig.perc & ncol(data) >= n.cutoff ){

    #print( paste( signature_name , "|" , "Weighted Mean" , sep=" " ) )

    gene <- intersect( rownames(data) , sig$gene_name)
    s <- sig[ sig$gene_name %in% gene, ]

    data <- data[ gene , ]

    geneSig <- NULL
    geneSig <- apply( data , 2 , function(x) ( exp(mean(x, na.rm=TRUE )) ) )
    names( geneSig ) <- colnames(data)

  }else{

    geneSig <- NA
    message("not enough samples and/or genes in a data")

  }

  return(geneSig)
}

#####################################################################
## Get signature score: PredictIO
#####################################################################

geneSigPredictIO <- function(dat.icb, clin = NULL, sig, sig.name, missing.perc, const.int =0.001, n.cutoff, sig.perc, study){

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

    geneSig <- NULL
    if( ncol(data)){

      #print( paste( signature_name , "|" , "Specific" , sep=" " ) )

      sensitive <- sig[ sig$weight == "sensitive" , ]$gene_name
      resistance <- sig[ sig$weight == "resistance", ]$gene_name

      IO_resistance <- NULL
      if( ifelse( is.null( nrow( data[ rownames(data) %in% resistance , ]) ) , 1 , nrow( data[ rownames(data) %in% resistance , ] ) ) / length( resistance ) > sig.perc ){
        gsvaPar <- gsvaParam(scale.fun( x=data ) , list(resistance))
        IO_resistance <- as.numeric(gsva(gsvaPar, verbose=FALSE))
      }

      IO_sensitive <- NULL
      if( ifelse( is.null( nrow( data[ rownames(data) %in% sensitive , ]) ) , 1 , nrow( data[ rownames(data) %in% sensitive , ] ) ) / length( sensitive ) > sig.perc ){
        gsvaPar <- gsvaParam(scale.fun( x=data ) , list(sensitive))
        IO_sensitive <- as.numeric(gsva(gsvaPar, verbose=FALSE))
      }

      if( !is.null( IO_resistance ) & !is.null( IO_sensitive ) ){

        geneSig <- IO_sensitive - IO_resistance
        names(geneSig) <- colnames(data)

          }else{

            geneSig <- NA
            message("not enough samples and/or genes in a data")

        }

    }else{

      geneSig <- NA
      message("not enough samples in a data")

    }

  return(geneSig)
}


#####################################################################
## Get signature score: COX_IS
#####################################################################

geneSigCOX_IS <- function(dat.icb, clin = NULL, sig, sig.name, missing.perc, const.int =0.001, n.cutoff, sig.perc, study){

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

    if( ifelse( is.null( nrow( data[ rownames(data) %in% sig$gene_name , ]) ) , 1 , nrow( data[ rownames(data) %in% sig$gene_name , ] ) ) / length( sig$gene_name ) > sig.perc & ncol(data) >= n.cutoff ){

      #print( paste( signature_name , "|" , "Specific" , sep=" " ) )

      geneSig <- NULL
      pos <- apply( data[ rownames(data) %in% sig[ sig$weight %in% 1 , ]$gene_name , ] , 2 , function(x){ ( sum( x) /  length( x) ) } )
      neg <- apply( data[ rownames(data) %in% sig[ sig$weight %in% -1 , ]$gene_name , ] , 2 , function(x){ ( sum( x) /  length( x) ) } )

      geneSig <- as.numeric( scale( pos / neg ) )
      names(geneSig) <- colnames(data)

      if( length( geneSig[ !is.nan( geneSig ) ] ) < n.cutoff ){

        geneSig <- NA
        message("not enough samples and/or genes in a data")

       }

    }else{

      geneSig <- NA
      message("not enough samples and/or genes in a data")

    }

  return(geneSig)
}

#####################################################################
## Get signature score: IPS
#####################################################################
## calculate Immunophenoscore

ipsmap <- function (x) {
  if (x<=0) {
    ips<-0
  } else {
    if (x>=3) {
      ips<-10
    } else {
      ips<-round(x*10/3, digits=0)
    }
  }
  return(ips)
}

IPS.fun <- function( expr, sig ){

  expr <- as.data.frame(expr)
  sample_names <- names(expr)

  unique_ips_genes <- as.vector(unique(sig$name))

  IPS <- NULL
  MHC <- NULL
  CP <- NULL
  EC <- NULL
  SC <- NULL
  AZ <- NULL

  # Gene names in expression file
  GVEC <- row.names( expr )
  # Genes names in IPS genes file
  VEC <- as.vector( sig$gene_name )
  # Match IPS genes with genes in expression file
  ind <- which( is.na( match( VEC , GVEC ) ) )

  if( length( ind ) ){
    sig <- sig[-ind,]
  }

  for (i in 1:length(sample_names)) {
    GE <- expr[[i]]
    mGE <- mean(GE, na.rm=TRUE)
    sGE <- sd(GE, na.rm=TRUE)
    Z1 <- (expr[as.vector(sig$gene_name),i] - mGE)/sGE
    W1 <- sig$weight
    coef <- NULL
    MIG <- NULL
    k <- 1

    for (gen in unique_ips_genes) {

      MIG[k] <- mean(Z1[which (as.vector(sig$name)==gen)], na.rm=TRUE)
      coef[k] <- mean(W1[which (as.vector(sig$name)==gen)], na.rm=TRUE)
      k <- k+1

    }

    WG <- MIG * coef
    MHC[i] <- mean(WG[1:10], na.rm=TRUE)
    CP[i] <- mean(WG[11:20], na.rm=TRUE)
    EC[i] <- mean(WG[21:24], na.rm=TRUE)
    SC[i] <- mean(WG[25:26], na.rm=TRUE)
    AZ[i] <- sum(MHC[i], CP[i], EC[i], SC[i])
    IPS[i] <- ipsmap(AZ[i])
  }

  AZ
}

geneSigIPS <- function(dat.icb, clin = NULL, sig, sig.name, missing.perc, const.int =0.001, n.cutoff, sig.perc, study){


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


    geneSig = NULL
    if( ncol(data) & nrow(data) > 10000 ){

        #print( paste( signature_name , "|" , "Specific" , sep=" " ) )

        geneSig <- as.numeric( scale( IPS.fun( expr = data , sig = sig) ) )
        names( geneSig ) <- colnames(data)

      }else{

      geneSig <- NA
      message("not enough samples and/or genes in a data")

      }

  return(geneSig)
}

#####################################################################
## Get signature score: IPRES
#####################################################################

geneSigIPRES <- function(dat.icb, clin = NULL, sig, sig.name, missing.perc, const.int =0.001, n.cutoff, sig.perc, study){

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

    IPRES.dat <- sig

    sig <- list()
    for( j in 1:length(IPRES.dat)){
      sig[[j]] <-  IPRES.dat[[j]]$gene_name
    }
    names(sig) <- names(IPRES.dat)

    #cancer_type <- names( table( dat_clin$cancer_type )[ table( dat_clin$cancer_type ) >= n.cutoff ] )
    #message(paste(study))
    #data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type ]
    data <- dat_expr
    remove <- rem(data, missing.perc, const.int)

    if( length(remove) ){
      data <- data[-remove,]
    }

      #print( paste( signature_name , "|" , "Specific" , sep=" " ) )
      scale.data <- scale.fun( data )
      geneSig <- NULL

      for(k in 1:length(sig)){

        if( ifelse( is.null( nrow( scale.data[ rownames(scale.data) %in% sig[[k]] , ]) ) , 1 , nrow( scale.data[ rownames(scale.data) %in% sig[[k]] , ] ) ) / length( sig[[k]] ) >= sig.perc & ncol(scale.data) >= n.cutoff ){

          geneSig[[k]] <- gsva(scale.data , list(sig[[k]]) , method = "ssgsea", verbose=FALSE)

         }else{

          geneSig[[k]] <- NA

         }

    }

    if( sum(is.na(geneSig)) > 0 & sum(is.na(geneSig)) != length(geneSig) ){ geneSig <- geneSig[!is.na(geneSig)]   }else{

        if(sum(is.na(geneSig)) > 0 & sum(is.na(geneSig)) == length(geneSig)){ geneSig <- NA }else{

          geneSig <- geneSig
        }
      }

      if( sum(!is.na(geneSig)) != 0){

        geneSig <- do.call(rbind, geneSig)
        geneSig <- apply( geneSig , 2 , mean , na.rm=TRUE )
        names(geneSig) <- colnames(data)

      }else{

        geneSig <- NA
        message("not enough samples and/or genes in a data")
      }

  return(geneSig)
}


#####################################################################
## Get signature score: PassON
#####################################################################

geneSigPassON <- function(dat.icb, clin = NULL, sig, sig.name, missing.perc, const.int =0.001, n.cutoff, sig.perc, study){

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

  PassON.dat <- sig

  sig <- list()
  for( j in 1:length(PassON.dat)){
    sig[[j]] <-  PassON.dat[[j]]$gene_name
  }
  names(sig) <- names(PassON.dat)

  #cancer_type <- names( table( dat_clin$cancer_type )[ table( dat_clin$cancer_type ) >= n.cutoff ] )
  #message(paste(study))
  #data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type ]
  data <- dat_expr
  remove <- rem(data, missing.perc, const.int)

  if( length(remove) ){
    data <- data[-remove,]
  }

  #print( paste( signature_name , "|" , "Specific" , sep=" " ) )
  scale.data <- scale.fun( data )
  geneSig <- NULL

  for(k in 1:length(sig)){

    if( ifelse( is.null( nrow( scale.data[ rownames(scale.data) %in% sig[[k]] , ]) ) , 1 , nrow( scale.data[ rownames(scale.data) %in% sig[[k]] , ] ) ) / length( sig[[k]] ) >= sig.perc & ncol(scale.data) >= n.cutoff ){

      geneSig[[k]] <- gsva(scale.data , list(sig[[k]]) , method = "ssgsea", verbose=FALSE)

    }else{

      geneSig[[k]] <- NA

    }

  }

  if( sum(is.na(geneSig)) > 0 & sum(is.na(geneSig)) != length(geneSig) ){ geneSig <- geneSig[!is.na(geneSig)]   }else{

    if(sum(is.na(geneSig)) > 0 & sum(is.na(geneSig)) == length(geneSig)){ geneSig <- NA }else{

      geneSig <- geneSig
    }
  }

  if( sum(!is.na(geneSig)) != 0){

    geneSig <- do.call(rbind, geneSig)
    geneSig <- apply( geneSig , 2 , mean , na.rm=TRUE )
    names(geneSig) <- colnames(data)

  }else{

    geneSig <- NA
    message("not enough samples and/or genes in a data")
  }

  return(geneSig)
}


#####################################################################
## Get signature score: IPSOV
#####################################################################

geneSigIPSOV <- function(dat.icb, clin = NULL, sig, sig.name, missing.perc, const.int =0.001, n.cutoff, sig.perc, study){

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

  IPSOV.dat <- sig

  sig <- list()
  for( j in 1:length(IPSOV.dat)){
    sig[[j]] <-  IPSOV.dat[[j]]$gene_name
  }
  names(sig) <- names(IPSOV.dat)

  #cancer_type <- names( table( dat_clin$cancer_type )[ table( dat_clin$cancer_type ) >= n.cutoff ] )
  #message(paste(study))
  #data <- dat_expr[ , dat_clin$cancer_type %in% cancer_type ]
  data <- dat_expr
  remove <- rem(data, missing.perc, const.int)

  if( length(remove) ){
    data <- data[-remove,]
  }

  #print( paste( signature_name , "|" , "Specific" , sep=" " ) )
  scale.data <- scale.fun( data )
  geneSig <- NULL

  for(k in 1:length(sig)){

    if( ifelse( is.null( nrow( scale.data[ rownames(scale.data) %in% sig[[k]] , ]) ) , 1 , nrow( scale.data[ rownames(scale.data) %in% sig[[k]] , ] ) ) / length( sig[[k]] ) >= sig.perc & ncol(scale.data) >= n.cutoff ){

      geneSig[[k]] <- gsva(scale.data , list(sig[[k]]) , method = "ssgsea", verbose=FALSE)

    }else{

      geneSig[[k]] <- NA

    }

  }

  if( sum(is.na(geneSig)) > 0 & sum(is.na(geneSig)) != length(geneSig) ){ geneSig <- geneSig[!is.na(geneSig)]   }else{

    if(sum(is.na(geneSig)) > 0 & sum(is.na(geneSig)) == length(geneSig)){ geneSig <- NA }else{

      geneSig <- geneSig
    }
  }

  if( sum(!is.na(geneSig)) != 0){

    geneSig <- do.call(rbind, geneSig)
    geneSig <- apply( geneSig , 2 , mean , na.rm=TRUE )
    names(geneSig) <- colnames(data)

  }else{

    geneSig <- NA
    message("not enough samples and/or genes in a data")
  }

  return(geneSig)
}
