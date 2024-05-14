## load libraries
library(survcomp)
library(survival)
library(GSVA)
library(data.table)
library(MultiAssayExperiment)

############################################################################
## Remove genes if with expression zero in 50% (or missing.perc) of sample
############################################################################

rem <- function(x, missing.perc, const.int){

  x <- as.matrix(x)
  x <- t(apply(x,1,as.numeric))

  # data is log2(TPM+0.001)
  r <- as.numeric(apply(x, 1, function(i) sum(round(i, 6) == round(log2(const.int), 6)) ))
  remove <- which(r > dim(x)[2]* missing.perc)
  return(remove)

 }

##################################################################################################
## Cox model: OS/PFS analyses and continuous expression, signature score, clinical data (var)
##################################################################################################

survCont <- function( surv , time , time.censor , var){

  data <- data.frame( surv=surv , time=time , variable=var )
  data <- data[!is.na(data$variable), ]
	data$time <- as.numeric(as.character(data$time))
	data$variable <- as.numeric( as.character(data$variable) )

  	for(i in 1:nrow(data)){

	    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){

	       data[ i , "time" ] <- time.censor
	       data[ i , "surv" ] <- 0

	       }
  	}

	cox <- coxph( Surv( time , surv ) ~ variable, data=data )

	hr <- summary(cox)$coefficients[, "coef"]
	se <- summary(cox)$coefficients[, "se(coef)"]
	n <- round(summary(cox)$n)
	low <- summary(cox)$conf.int[, "lower .95"]
	up <- summary(cox)$conf.int[, "upper .95"]
	pval <- summary(cox)$coefficients[, "Pr(>|z|)"]

	res <- c( hr , se , n, low , up , pval )
	names(res) <- c( "HR" , "SE" , "N", "Low" , "Up" , "Pval")

	return(res)

 }

########################################################################################
## Cox model: OS/PFS analyses and binary expression or signature score (low vs high)
########################################################################################
# n0.cutoff: minimum number of samples less than cutoff
# n1.cutoff: minimum number of samples greater than cutoff
# var.type: if variable (var) is dicho (by default), then var.type is TRUE

survDicho <- function(surv , time , time.censor , var , n0.cutoff, n1.cutoff, method ="median", var.type = TRUE){

  data <- data.frame( surv=surv , time=time , variable=var )
  data <- data[!is.na(data$variable), ]
  data$time <- as.numeric(as.character(data$time))

  if(var.type != TRUE){

    if( method == "median"){
      data$variable <- ifelse( as.numeric(as.character(data$variable)) >= median(as.numeric(as.character(data$variable))) , 1 , 0 )
    }

    if( method == "Q1" ){
      data$variable <- ifelse( as.numeric(as.character(data$variable)) >= quantile(as.numeric(as.character(data$variable)))["25%"] , 1 , 0 )
    }

    if( method == "Q3" ){
      data$variable <- ifelse( as.numeric(as.character(data$variable)) >= quantile(as.numeric(as.character(data$variable)))["75%"] , 1 , 0 )
    }

  }

  for(i in 1:nrow(data)){

    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){
      data[ i , "time" ] = time.censor
      data[ i , "surv" ] = 0

    }
  }

  if( length( data$variable[ data$variable == 1 ] )>= n1.cutoff &
      length( data$variable[ data$variable == 0 ] ) >= n0.cutoff ){

    cox <- coxph( formula= Surv( time , surv ) ~ variable , data=data )
    hr <- summary(cox)$coefficients[, "coef"]
    se <- summary(cox)$coefficients[, "se(coef)"]
    n <- round(summary(cox)$n)
    low <- summary(cox)$conf.int[, "lower .95"]
    up <- summary(cox)$conf.int[, "upper .95"]
    pval <- summary(cox)$coefficients[, "Pr(>|z|)"]

  } else{

    hr <- NA
    se = NA
    low <- NA
    up <- NA
    pval <- NA

  }

  res <- c( hr , se , n, low , up , pval )
  names(res) <- c( "HR" , "SE" , "N", "Low" , "Up" , "Pval")
  return(res)

 }

##################################################################################################
## Kaplan-Meier (KM) plot: OS/PFS analyses and binary expression or signature score (low vs high)
##################################################################################################
# n0.cutoff: minimum number of samples less than cutoff
# n1.cutoff: minimum number of samples greater than cutoff
# var.type: if variable (var) is dicho (by default), then var.type is TRUE

KMPlot <- function( surv , time , time.censor , var , title , xlab, ylab,
                    method = "median", n0.cutoff, n1.cutoff, var.type = TRUE){

  data <- data.frame( surv=surv , time=time , variable=var )
  data <- data[!is.na(data$variable), ]
  data$time <- as.numeric(as.character(data$time))

  if( var.type != TRUE){

    if( method == "median"){
      bin.cutoff <- median(as.numeric(as.character(data$variable)))
      data$variable <- ifelse( as.numeric(as.character(data$variable)) >= bin.cutoff , 1 , 0 )
    }

    if( method == "Q1" ){
      bin.cutoff <- quantile(as.numeric(as.character(data$variable)))["25%"]
      data$variable <- ifelse( as.numeric(as.character(data$variable)) >= bin.cutoff , 1 , 0 )
    }

    if( method == "Q3" ){
      bin.cutoff <- quantile(as.numeric(as.character(data$variable)))["75%"]
      data$variable <- ifelse( as.numeric(as.character(data$variable)) >= bin.cutoff , 1 , 0 )
    }

  }

  for(i in 1:nrow(data)){

    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){
      data[ i , "time" ] = time.censor
      data[ i , "surv" ] = 0

    }
  }

  if( length( data$variable[ data$variable == 1 ] )>= n1.cutoff & length( data$variable[ data$variable == 0 ] ) >= n0.cutoff ){

    km.coxph.plot( Surv( time, surv ) ~ variable , data.s = data, x.label = xlab, y.label = ylab,
                   main.title = paste( title , "\n(cutoff=" , round( bin.cutoff , 2 ) , ")" , sep="" ) ,
                   sub.title = "",
                   leg.text = c( "Low" , "High"),
                   leg.pos = "topright",
                   .col = c( "#b2182b","#2166ac"),
                   show.n.risk = TRUE,
                   n.risk.step = 7,
                   n.risk.cex = 0.70,
                   ylim = c(0,1),
                   leg.inset = 0,
                   .lwd = 2 ,
                   verbose=FALSE )
   }

}


##########################################################################################
## volcano plot for signatures or genes association with immunotherapy responses results
##########################################################################################

volcanoPlot <- function(feature, coef, pval, padj, pos.cutoff, neg.cutoff, x.lab, padj.label, cutoff){

  data <- data.frame(feature = feature,
                     coef = coef,
                     pval = pval,
                     FDR = padj)

  if( padj.label == FALSE){

    data$diffexpressed <- "NO"
    data$diffexpressed[data$coef > 0 & data$pval < cutoff] <- paste(paste("Pval < ", cutoff, sep=""), "Coef > 0", sep=", ")
    data$diffexpressed[data$coef < 0 & data$pval < cutoff] <- paste(paste("Pval < ", cutoff, sep=""), "Coef < 0", sep=", ")

    mycolors <- c( "#984ea3","#386cb0", "#999999")
    names(mycolors) <- c(paste(paste("Pval < ", cutoff, sep=""), "Coef > 0", sep=", "),
                         paste(paste("Pval < ", cutoff, sep=""), "Coef < 0", sep=", "),
                         "NO")

    data$delabel <- NA
    data <- data[order(data$pval, decreasing = FALSE), ]
    id_pos <- data[data$coef > 0 , "feature"][1:pos.cutoff]
    id_neg <- data[data$coef < 0 , "feature"][1:neg.cutoff]
    id <- c(id_pos, id_neg)

    for(j in 1:length(id)){
      k <- which(data$feature == id[j])
      data$delabel[k] <- data[k, ]$feature
    }

  }else{

    data$diffexpressed <- "NO"
    data$diffexpressed[data$coef > 0 & data$FDR < cutoff] <- paste(paste("FDR < ", cutoff, sep=""), "Coef > 0", sep=", ")
    data$diffexpressed[data$coef < 0 & data$FDR < cutoff] <- paste(paste("FDR < ", cutoff, sep=""), "Coef < 0", sep=", ")

    mycolors <- c( "#984ea3","#386cb0", "#999999")
    names(mycolors) <- c(paste(paste("FDR < ", cutoff, sep=""), "Coef > 0", sep=", "),
                         paste(paste("FDR < ", cutoff, sep=""), "Coef < 0", sep=", "),
                         "NO")

    data$delabel <- NA
    data <- data[order(data$FDR, decreasing = FALSE), ]
    id_pos <- data[data$coef > 0 , "feature"][1:pos.cutoff]
    id_neg <- data[data$coef < 0 , "feature"][1:neg.cutoff]
    id <- c(id_pos, id_neg)

    for(j in 1:length(id)){
      k <- which(data$feature == id[j])
      data$delabel[k] <- data[k, ]$feature
    }

  }

  ggplot(data=data, aes(x=coef, y=-log10(pval), col= diffexpressed)) +
    geom_point(size = 3) + theme_minimal() +
    ylab("-log10 P value") +
    xlab(x.lab) +
    scale_colour_manual(values = mycolors) +
    theme(
      axis.text.x=element_text(size=12,  face="bold"),
      axis.title=element_text(size=12,face="bold"),
      axis.text.y=element_text(size=12, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position="bottom",
      legend.text = element_text(size = 9, face="bold"),
      legend.title = element_blank()) +
      geom_text_repel(aes(label= delabel),
                    size = 2.7,
                    color = "black",
                    min.segment.length = 0,
                    na.rm = TRUE, direction = "both",
                    seed = 2356,
                    fontface= "bold",
                    max.overlaps = 50)
  }

