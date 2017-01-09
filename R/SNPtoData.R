#' Generate Omics simulated data
#'
#' Generate Omics simulated data based on Glycan genetics and Measurement error.
#'
#' @author
#' MIMOmics consortium
#'
#' Maintainer
#' MIMOmics (\email{mimomics@lumc.nl})
#'
#' @section Background:
#' One.
#'
#' @section Overview:
#' Two.
#'
#' @section Details:
#' Three.
#'
#' @docType package
#' @name MIMOmicsData
#' @keywords MIMOmics Omics DataPackage
#' @importFrom MASS mvrnorm
#' @importFrom stats complete.cases quantile rbinom rnorm runif
#' @importFrom utils read.csv read.table
#' @import dplyr gplots stringr magrittr
NULL

# library("rex") # for one_of
# library(dplyr) # for data.frame manipulation
# library(gplots) # for heatmap.2
# library(stringr) # for strings
# library(magrittr) # for the %>% operator
# also Need MASS::mvrnorm

#' Calculate Sum of Squares
#'
#' @param X Numeric vector or matrix.
#' @return The sum of squared elements of \eqn{X}
#' @details This is the Frobenius norm of \eqn{X}.
#' @examples
#' ssq(tcrossprod(1:5))
#' ssq(rnorm(1e5))/1e5
#' @export
ssq <- function(X) {
  return(sum(X^2))
}


#' Escoufier's RV coeff, measure of similarity
#'
#' @param X Numeric Matrix
#' @param Y Numeric Matrix
#' @details This is an implementation of Escoufier's RV coefficient.
#' The memory needed is of order \eqn{p^2} with \eqn{p} the number of variables.
#' @export
coeffRV <- function(X,Y){
  X=scale(X,scale=F)
  Y=scale(Y,scale=F)
  ssq(t(Y)%*%X) / sqrt(ssq(t(X)%*%X)*ssq(t(Y)%*%Y))
}

# ############## SET PARAMETERS
# #######################################################################################
# #Number of X2 variables
# p_X2 = 100
# # Set a proportion between 0 and 1 (corresponding to RV coeff between 0 and 0.98, depends on the weight_shared_SNPs)
# # This is most influential on the cor between X1 and X2, determines to which extent variables in X1 and X2 agree
# similarity_data = 1
# # Set within correlation in X2 by changing proportion of SNPs used in each X2 variable.
# within_corr = 0.1
# # Set weight of the shared SNPs for X2. More weight means shared SNPs are more important for X2
# #       However as the linear combinations in X1 are not the same as the linear combinations in X2,
# #       increasing the weights does not necessarily increase the similarity
# #       Set between 0 and Inf, 0 gives approx independent X1 and X2 if similarity_data = 0
# weight_shared_SNPs = 0
# # Set Measurement error correlations and variances for X1 (first element cor, then var)
# meas_add_X1 = c(0.25, 0.5)
# meas_mult_X1 = c(0.05, 0.1)
# meas_add_X2 = c(0.5, 1) / 4
# meas_mult_X2 = c(0.01, 0.02)
# #######################################################################################

#' Generate simulated omics data
#'
#' Generate simulated omics data using genetics from the MIMOmics consortium.
#'
#' @param nr_X2vars Number of X2 variables
#' @param similarity_data Between 0 and 1. How similar are X1 and X2?
#' @param within_corr Between 0 and 1. Within correlation of X2.
#' @param weight_shared_SNPs Larger or equal to zero. Obsolete!
#' @param meas_add_X1 A vector of two numerics.
#' @param meas_mult_X1 A vector of two numerics.
#' @param meas_add_X2 A vector of two numerics.
#' @param meas_mult_X2 A vector of two numerics.
#'
#' @return A list with the data:
#' \item{X1}{The first dataset without measurement error.}
#' \item{X1m}{The first dataset with measurement error.}
#' \item{X2}{The second dataset without measurement error.}
#' \item{X2m}{The first dataset with measurement error.}
#' \item{y}{Continuous outcome}
#' \item{d1}{Binary outcome ...}
#' \item{d}{Binary outcome ...}
#'
#' @details Number of X2 variables
#' p_X2 = 100
#' Set a proportion between 0 and 1 (corresponding to RV coeff between 0 and 0.98, depends on the weight_shared_SNPs)
#' This is most influential on the cor between X1 and X2, determines to which extent variables in X1 and X2 agree
#' similarity_data = 1
#' # Set within correlation in X2 by changing proportion of SNPs used in each X2 variable.
#' within_corr = 0.1
#' # Set weight of the shared SNPs for X2. More weight means shared SNPs are more important for X2
#' #       However as the linear combinations in X1 are not the same as the linear combinations in X2,
#' #       increasing the weights does not necessarily increase the similarity
#' #       Set between 0 and Inf, 0 gives approx independent X1 and X2 if similarity_data = 0
#' weight_shared_SNPs = 0
#' # Set Measurement error correlations and variances for X1 (first element cor, then var)
#' meas_add_X1 = c(0.25, 0.5)
#' meas_mult_X1 = c(0.05, 0.1)
#' meas_add_X2 = c(0.5, 1) / 4
#' meas_mult_X2 = c(0.01, 0.02)
#' @export

GenData <- function(nr_X2vars = 100,
                    similarity_data = 1,
                    within_corr = 0.1,
                    weight_shared_SNPs = 0,
                    meas_add_X1 = c(.25, .5),
                    meas_mult_X1 = meas_add_X1/5,
                    meas_add_X2 = meas_add_X1/2,
                    meas_mult_X2 = meas_mult_X1/5
                    )
{
  p_X2 = nr_X2vars
  SNPkor = read.table(file = system.file("extdata","data_Kor.raw",package='MIMOmicsData'),header=TRUE)
  SNPvis = read.table(file = system.file("extdata","data_Vis.raw",package='MIMOmicsData'),header=TRUE)

  ## SELECT intersecting names
  namesBoth = intersect(names(SNPkor),names(SNPvis))
  SNPkor %>% select(one_of(namesBoth)) -> SNPkor2
  SNPvis %>% select(one_of(namesBoth)) -> SNPvis2

  # Complete cases
  SNPkor2 %>% filter(complete.cases(.)) -> SNPkor3
  SNPvis2 %>% filter(complete.cases(.)) -> SNPvis3

  #cat(paste(as.character(SNPvis$FID[unique(which(is.na(SNPvis2),arr.ind = T)[,1])]),as.character(SNPvis$FID[unique(which(is.na(SNPvis2),arr.ind = T)[,1])])),sep='\n')

  ## REMOVE Heterozygous variables
  SNPkor3 %>% select(-c(1,3:6)) %>% rename(ID = IID) %>% select(-ends_with("HET")) -> SNPkor4
  SNPvis3 %>% select(-c(1,3:6)) %>% rename(ID = IID) %>% select(-ends_with("HET")) -> SNPvis4

  ## RBIND the data
  SNPkor4 %>% rbind(SNPvis4) -> SNPdat

  # SNPp = list met p-values
  hits_glycans = as.list(read.csv(system.file("extdata","UPLC_IgG_varianceExplained_byGlycan.csv",package='MIMOmicsData'), header=T, sep=';'))
  SNPhits = sapply(1:length(hits_glycans[[2]]), function(i) strsplit(as.character(hits_glycans$SNPs[i]),split=';'))
  SNPp = sapply(1:length(hits_glycans[[2]]), function(i) hits_glycans$beta[i] %>% as.character %>% strsplit(split=';') %>% extract2(1) %>% as.numeric)

  SNPhits2 = sapply(2:ncol(SNPdat), function(k) sapply(1:length(SNPhits), function(i) max(grep(str_sub(names(SNPdat)[k],end=-3), str_to_upper(str_trim(SNPhits[[i]]))),0)))
  colnames(SNPhits2) = colnames(SNPdat)[-1]

  # WeightMat1 contains the weights for X1
  SNPscore = matrix(NA, nrow(SNPdat), nrow(SNPhits2))
  WeightMat1 = matrix(0, ncol(SNPdat[,-1]), nrow(SNPhits2))
  for(IGP in 1:nrow(SNPhits2)){
    ind = which(SNPhits2[IGP,] > 0)
    WeightMat1[ind,IGP] = SNPp[[IGP]][SNPhits2[IGP,ind]]
    SNPscore[,IGP] = (as.matrix(SNPdat[,names(ind)]) %*% SNPp[[IGP]][SNPhits2[IGP,ind]])
  }
  SNPscore %<>% as.data.frame
  names(SNPscore) <- hits_glycans$glycan
  X1 = as.matrix(SNPscore)
  colnames(X1) <- names(SNPscore)
  rownames(X1) <- rownames(SNPscore)

  # data2b contain random SNPs
  data2b <- read.table(file = system.file("extdata","datareduce.raw",package='MIMOmicsData'),header=TRUE)[,-(2:6)]
  # order so that names match with first dataset SNPdat
  ordernames2b <- order(sapply(1:nrow(data2b), function(i) which(as.character(data2b$FID[i]) == as.character(SNPdat$ID))))
  data2b<- data2b[ordernames2b,-1]

  # data2 contains both first set of SNPs and random SNPs
  data2 = cbind(SNPdat, data2b)
  Pval = sort(unlist(lapply(SNPp,function(e) e)))
  # Give first SNPs higher weights, sampled from "observed" p-values
  WeightMat = weight_shared_SNPs*matrix(sample(Pval[1:200],29*p_X2,replace = T)*rbinom(29*p_X2,1,within_corr),29,p_X2)
  WeightMat = rbind(WeightMat, matrix(sample(Pval[-(1:200)],147*p_X2,replace = T)*rbinom(147*p_X2,1,within_corr),147,p_X2))
  # This X2-variable depends on all SNPs
  very_important_X2var = sample(1:ncol(WeightMat),1)
  WeightMat[,very_important_X2var] = sample(Pval,147+29,replace = T)
  WeightMat = WeightMat / sqrt(ssq(WeightMat)) * 100
  # Similarity as convex combi of independent and same SNPs
  permuted_cols = sample(1:72)
  WeightMat[1:29,1:72] = (1-similarity_data) * WeightMat[1:29,1:72]*0 + similarity_data * WeightMat1[,permuted_cols]/ sqrt(ssq(WeightMat1)) * 100
  WeightMat[1:29,73:p_X2] %<>% multiply_by(0)
  X2 = as.matrix(data2[,-1]) %*% WeightMat

  # Meas.err.
  covUa <- matrix(meas_add_X1[1], ncol(SNPscore),ncol(SNPscore)) # additive
  diag(covUa) = meas_add_X1[2]
  covUm <- matrix(meas_mult_X1[1], ncol(SNPscore),ncol(SNPscore)) # multiplicative
  diag(covUm) = meas_mult_X1[2]
  Ua <- MASS::mvrnorm(nrow(SNPscore), rep(0, nrow(covUa)), covUa)
  Um <- MASS::mvrnorm(nrow(SNPscore), rep(0, nrow(covUa)), covUm)
  X1m <- X1/sqrt(ssq(X1)/ssq(exp(Um) + Ua)/10) * exp(Um) + Ua

  colnames(X2) <- paste("V",1:ncol(X2),sep="")
  rownames(X2) <- rownames(SNPscore)
  covUa <- matrix(meas_add_X2[1], ncol(X2),ncol(X2)) # additive
  diag(covUa) = meas_add_X2[2]
  covUm <- matrix(meas_mult_X2[1], ncol(X2),ncol(X2)) # multiplicative
  diag(covUm) = meas_mult_X2[2]
  Ua <- MASS::mvrnorm(nrow(X2), rep(0, nrow(covUa)), covUa)
  Um <- MASS::mvrnorm(nrow(X2), rep(0, nrow(covUa)), covUm)
  X2m <- X2/sqrt(ssq(X2)/ssq(exp(Um) + Ua)/10) * exp(Um) + Ua

  #Add outcomes#

  #Association with outcome via a subset of the 29 snps which generate X1 (not added value of X2)#
  p1<-15
  beta1<-c(runif(p1,0,1),rep(0,29-p1))
  y1<-as.matrix(SNPdat[,-1])%*%beta1+rnorm(nrow(SNPdat),0,1) #continuous
  d1<-(y1>=quantile(y1)[4])                                  #binary

  #Association with outcome via a subset of the 149 snps which generate X2 (and not X1)#
  p2<-15
  beta2<-c(runif(p2,0,1),rep(0,147-p2))
  y2<-as.matrix(data2b)%*%beta2+rnorm(nrow(data2b),0,1)

  #Association with outcome via both X1 and X2 (added value of X2)#
  y<-y1+y2
  d<-(y1>=quantile(y1)[4])

  outp = list(X1 = X1, X1m = X1m, X2 = X2, X2m = X2m, y = y, d1 = d1, d = d)
  class(outp) <- "MIMOmicsData"
  return(outp)
}

#' Plot heatmaps of simulated data
#'
#' @param x A MIMOmicsData object as generated by \code{\link{GenData}}.
#' @return NULL
#' @keywords internal
#' @export
plot.MIMOmicsData <- function(x, ...){
  with(x, {
    suppressWarnings(heatmap.2(cor(X1[,-which(colSums(X1)==0)]), Rowv=F, Colv=F, col = bluered, breaks=seq(-1,1,length.out = 100),symm=T, trace="n",main='X1 corr w/o zero-sd IGPs'))
    suppressWarnings(heatmap.2(cor(X1m), Rowv=F, Colv=F, col = bluered, breaks=seq(-1,1,length.out = 100),symm=T, trace="n",main='X1m with all IGP'))
    suppressWarnings(heatmap.2(cor(X2), Rowv=F, Colv=F, col = bluered, breaks=seq(-1,1,length.out = 100),symm=T, trace="n",main='X2'))
    suppressWarnings(heatmap.2(cor(X2m), Rowv=F, Colv=F, col = bluered, breaks=seq(-1,1,length.out = 100),symm=T, trace="n",main='X2m'))
    if(ncol(X2m)>=100){
      suppressWarnings(heatmap.2(cor(X1m,X2m[,1:100]), Rowv=F, Colv=F, col = bluered, breaks=seq(-1,1,length.out = 100),symm=F, trace="n",main='Cor between X1m and X2m (first 100 variables)'))
    }
  })
}

#' Print RV coefs of simulated data
#'
#' @inheritParams plot.MIMOmicsData
#' @return NULL
#' @keywords internal
#' @export
print.MIMOmicsData <- function(x, digits = 3, ...){
  with(x, {
    cat("RV coeff between error free X1 and X2:", coeffRV(X1, X2) %>% round(digits),'\n')
    cat("RV coeff between noisy X1m and X2m:", coeffRV(X1m, X2m) %>% round(digits),'\n')
    cat("RV coeff between X1 and X1m:", coeffRV(X1, X1m) %>% round(digits),'\n')
    cat("RV coeff between X2 and X2m:", coeffRV(X2, X2m) %>% round(digits),'\n')
  })
}
