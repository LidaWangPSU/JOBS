#' Title
#' @param eqtls
#' @param eqtls_se
#' @param weights
#' @return
#' @export
#'
#' @examples jobs.function(eqtls,eqtls_se,weight)
jobs.eqtls<-function(beta,se,weight){
  df<-beta
  df_s<-se

  cat(paste0("Check NA in sc-eqtl \n"))
  na.list<-list()
  for(r in 3:ncol(df)){
    snps<-df[which(is.na(df[,r])),1]
    names<-colnames(df)[r]
    na.list[[names]]<-snps
    cat(paste0(length(snps)," NA in ",names," \n"))
  }

  df[is.na(df)]<-0
  df_s[is.na(df_s)]<-1

  x_old<-df[,3:ncol(df)]

  cat(paste0("Check if all sc-eqtl effects are 0 \n"))
  cat(paste0(nrow(x_old), " gene-snp pairs in total \n"))

  sum1<-rowSums(x_old)

  length(sum1)
  length(which(sum1!=0))

  df<-df[which(sum1!=0),]

  df_s<-df_s[which(sum1!=0),]

  cat(paste0(nrow(df), " gene-snp pairs passed the filtering \n"))

  weight<-t(as.matrix(weight))
  weight_total<-t(matrix(as.numeric(weight),nrow=ncol(x_old),ncol=nrow(x_old)))
  pred<-as.matrix(x_old)%*%t(weight)
  sum1<-rowSums(x_old)
  se<-df_s[,3:ncol(df)]
  se_new<-se
  x_new<-x_old


  cat(paste0("Step2: Start analyzing ",nrow(x_old)," gene-snp pairs \n"))
  ###jointly model eqtls
  for(j in 1:ncol(x_old)){
    #print(j)
    ###estimate se by fisher information
    se_new[,j]<-sqrt(df_s[,2]^2*se[,j]^2/(df_s[,2]^2+weight_total[,j]^2*se[,j]^2))
    ###estimate effect size by maximize likelihood jointly
    x_new[,j]<-(x_old[,j]*df_s[,2]^2+df[,2]*weight_total[,j]*se[,j]^2-(pred-x_old[,j]*weight_total[,j])*weight_total[,j]*se[,j]^2)/(df_s[,2]^2+se[,j]^2*weight_total[,j]^2)
  }

  ##check results
  head(x_new)
  head(x_old)

  head(se_new)
  head(se)

  df_new<-as.matrix(cbind(df[,1:2],x_new))
  df_s_new<-as.matrix(cbind(df_s[,1:2],se_new))


  colnames(df_new)<-colnames(df)
  colnames(df_s_new)<-colnames(df_s)

  ###replace NA list with NA
  for(r in 3:ncol(df_s)){
    names<-colnames(df)[r]
    snps<-na.list[[names]]
    df_new[na.omit(match(snps,df_new)),r]<-NA
    df_s_new[na.omit(match(snps,df_new)),r]<-NA
  }

  cat("Finished! \n")

  res<-list("eqtls_new"=as.data.frame(df_new),"eqtls_se_new"=as.data.frame(df_s_new))

  return(res)
}


#' Title
#' @param eqtls
#' @param eqtls_se
#' @return
#' @export
#'
#' @examples jobs.function(eqtls,eqtls_se,weights)
jobs.nnls.weights<-function(beta,se){
  df<-as.data.frame(beta)
  df_s<-as.data.frame(se)
  ###check NA
  #cat(paste0("Check NA in sc-eqtl \n"))
  na.list<-list()
  for(r in 3:ncol(df_s)){
    snps<-df[which(is.na(df[,r])),1]
    names<-colnames(df)[r]
    na.list[[names]]<-snps
    #cat(paste0(length(snps)," NA in ",names," \n"))
  }

  df[is.na(df)]<-0
  df_s[is.na(df_s)]<-1
  ###check if all sc-eqtl effects are 0
  x_old<-df[,3:ncol(df)]
  #cat(paste0("Check if all sc-eqtl effects are 0 \n"))
  #cat(paste0(nrow(x_old), " gene-snp pairs in total \n"))
  sum1<-rowSums(abs(x_old))
  df<-df[which(sum1!=0),]
  #cat(paste0(nrow(df), " gene-snp pairs passed the filtering \n"))
  df_s<-df_s[which(sum1!=0),]

  ########estimate frac by NNLS
  pair<-df[,1]
  Y<-df[,2]
  X<-as.matrix(df[,3:ncol(df)])

  cat(paste0("Step1: Estimating cell type weights using NNLS \n"))
  library(pracma)
  Aeq <- matrix(rep(1, ncol(X)), nrow= 1)
  beq <- c(1)

  # Lower and upper bounds of the parameters, i.e [0, 1]
  lb <- rep(0, ncol(X))
  ub <- rep(1, ncol(X))

  # And solve:
  frac<-lsqlincon(X, Y, Aeq= Aeq, beq= beq, lb= lb, ub= ub)


  frac<-t(as.data.frame(frac))

  colnames(frac)<-colnames(X)

  weight<-frac

  weight<-as.numeric(weight)

  return(weight)
}


