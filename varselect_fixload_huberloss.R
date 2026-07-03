

FLCM.select.huber<-function(y,mydata,nbasis1,nbasis2,phiseq,alphas,lambdas)
{
GrandY<-y
n=length(unique(mydata$id))
T<-list()
for(i in 1:n)
{ind<-which(mydata$id==i)
T[[i]]<-mydata$time[ind]
}
name<-names(mydata)
name<-name[-c(1,2)]
p<- ncol(mydata)-2
UT<-c()
for( i in 1:n)
{UT<-union(UT,T[[i]])}
UT<-sort(UT)
knots = seq(min(UT),max(UT),l=(nbasis2-2))
norder = 4
nbasis = length(knots)- 2+ norder 
dayrng = c(min(UT),max(UT))

bbasis = create.bspline.basis(dayrng,nbasis,norder,knots)

bbasisMat<-list()
for (i in 1:n)
{bbasisMat[[i]] = eval.basis(T[[i]],bbasis)}

in.mat = inprod(bbasis,bbasis) #penalty matrix for function，计算bbasis与bbasis的内积矩阵

in.derivmat=inprod(bbasis,bbasis,int2Lfd(2),int2Lfd(2)) #penalty matrix for 2nd deriv.

basismat<-bbasisMat
for (i in 1:n){
  colnames(basismat[[i]])<-NULL}
XMatstar<-list()
for(i in 1:n)
{XMatstar[[i]]<-array(0,dim=c(p,length(T[[i]]),nbasis))
for (j in 1:p){
  for(l in 1:length(T[[i]])){
    indtemp<-which(mydata$id==i)
    datatemp<-mydata[indtemp,]
    datatempj<-datatemp[,2+j]
    XMatstar[[i]][j,l,]=datatempj[l]*basismat[[i]][l,]
   
  }
}
}
knotsint = seq(min(UT),max(UT),l=(nbasis1-2))
norderint = 4
nbasisint = length(knotsint) + norderint - 2
dayrngint = c(min(UT),max(UT))
bbasisint = create.bspline.basis(dayrngint,nbasisint,norderint,knotsint)

bbasisMatint<-list()
for (i in 1:n)
{bbasisMatint[[i]] = eval.basis(T[[i]],bbasisint)
}
B<-bbasisMatint
GrandB<-NULL
for(i in 1:n)
{GrandB<-rbind(GrandB,B[[i]])

}
group<-c()
d<-nbasis-nbasisint
for (j in 1:(p+1))
{
  if(j==1)
  {group[((nbasisint*(j-1))+1):(nbasisint*j)]<-0}

  if(j==2)
  {group[((nbasis*(j-1))+1-d):((nbasis*j)-d)]<-0}
  if(j>2)
  {group[((nbasis*(j-1))+1-d):((nbasis*j)-d)]<-j-2}#
  #################
}

ZMatstar<-list()
for(i in 1:n)
{ZMatstar[[i]]<-array(0,dim=c(p,length(T[[i]]),nbasis))
for (j in 1:p){
  for(l in 1:length(T[[i]])){
    ZMatstar[[i]][j,l,]=XMatstar[[i]][j,l,]
  }
}
}
Zimatlist<-list()
for (i in 1:n)
{Zimatlist[[i]]<-matrix(0,length(T[[i]]),nbasis*p)
for(l in 1:length(T[[i]]))
{for(j in 1 :p){
  Zimatlist[[i]][l,((nbasis*(j-1))+1):(nbasis*j)]<-ZMatstar[[i]][j,l,]
 
}
}
}
GrandZ<-NULL
for(i in 1:n)
{GrandZ<-rbind(GrandZ, Zimatlist[[i]])
}
group_no_intercept<-group
group<-c(0,group)
beta <- rep(1, length(group))
beta_new <- rep(2, length(group))
beta_1 <- beta
delta <- 1.5
epsilon <- 10^(-6)
max_iter <- 10000
iter_operation<-function(X, y, beta, delta,epsilon,beta_new,beta_1,group,alpha,lambda){
  k<-1
  while (norm(beta_new - beta_1, "2") > epsilon)
  {
    residuals <- y - X %*% beta
    gradient <- -t(X) %*% ifelse(abs(residuals) <= delta, residuals, delta * sign(residuals)) / nrow(X)
    
    beta_temp <- beta - alpha * gradient
    for (i in unique(group)[-1]) {
      group_indices<-which(group==i)
      norm_g <- norm(beta_temp[group_indices, , drop = FALSE], "2")
      beta_temp[group_indices, ] <- pmax(1 - (lambda * alpha) / norm_g, 0) * beta_temp[group_indices, ]

    }
    
    beta_new <- beta_temp
    
    beta_1<-beta
    beta <- beta_new
    k<-k+1
    if (k>max_iter){
      break
    }
  }
  print(paste("alpha=",alpha," lambda=",lambda))
  return(beta)
}



library(parallel)
ncores <-detectCores()-1
cl <- makeCluster(ncores)
cvphi<-function(phi){

  Kphi<-in.mat+phi*in.derivmat
  R<-chol(Kphi)
  L<-t(R)
  M<-solve(R)
  ZMatstar<-list()
  for(i in 1:n)
  {ZMatstar[[i]]<-array(0,dim=c(p,length(T[[i]]),nbasis))
  for (j in 1:p){
    for(l in 1:length(T[[i]])){
      ZMatstar[[i]][j,l,]=XMatstar[[i]][j,l,]%*%M

    }
  }
  }
  Zimatlist<-list()#与上一个Zimatstar计算流程一致
  for (i in 1:n)
  {Zimatlist[[i]]<-matrix(0,length(T[[i]]),nbasis*p)
  for(l in 1:length(T[[i]]))
  {for(j in 1 :p){
    Zimatlist[[i]][l,((nbasis*(j-1))+1):(nbasis*j)]<-ZMatstar[[i]][j,l,]

  }
  }
  }
  GrandZ<-NULL

  for(i in 1:n)
  {GrandZ<-rbind(GrandZ, Zimatlist[[i]])
  }
  Grandzz<-cbind(rep(1,nrow(GrandB)),GrandB,GrandZ)#第一列用来作为截距

  X<-Grandzz


  
  perform_metrics<- array(NA, dim = c(length(alphas), length(lambdas),6))
  
  gammlist<- array(NA, dim = c(length(alphas), length(lambdas),length(group)))
  y_pre_list<- array(NA, dim = c(length(alphas), length(lambdas),length(y)))
  for (a in 1:length(alphas)) {
      for (l in 1:length(lambdas)) {
          beta_train <- iter_operation(X, y, beta, delta,epsilon,beta_new,beta_1,group,alphas[a],lambdas[l])
          y_train<-X %*% beta_train
          # 计算性能指标，这里是Ebic准则
          ny<-length(y_train)
          nk<-length(which(beta_train!=0))
          RSS<-rmse(y,y_train)^2*ny
          deltahat2<-RSS/(ny-nk)
          lnL<-(-ny/2)*log(2*pi)-(ny/2)*log(deltahat2)-RSS/(2*deltahat2)
          EBIC<--2*lnL+nk*log(ny)+nk*(nk+1)/(2*ny)
          BIC<--2*lnL+nk*log(ny)
          AIC<--2*lnL+2*nk
          
          perform_metrics[a,l,]<-c(EBIC,BIC,AIC,mape(y,y_train),mae(y,y_train),rmse(y,y_train))
          #6个表现指标，都是越小越好
          gammlist[a, l,]<-beta_train
          y_pre_list[a, l,]<-y_train

      }
  }
  
  best_ebic_index<-which(perform_metrics[,,1] == min(perform_metrics[,,1]), arr.ind = TRUE)
  best_bic_index<-which(perform_metrics[,,2] == min(perform_metrics[,,2]), arr.ind = TRUE)
  best_aic_index<-which(perform_metrics[,,3] == min(perform_metrics[,,3]), arr.ind = TRUE)
  best_mape_index<-which(perform_metrics[,,4] == min(perform_metrics[,,4]), arr.ind = TRUE)
  best_mae_index<-which(perform_metrics[,,5] == min(perform_metrics[,,5]), arr.ind = TRUE)
  best_rmse_index<-which(perform_metrics[,,6] == min(perform_metrics[,,6]), arr.ind = TRUE)
  
  # 找到最佳的alpha和lambda
  best_alpha_ebic <- alphas[best_ebic_index[1, 1]]
  best_lambda_ebic <- lambdas[best_ebic_index[1, 2]]
  # 最终模型
  gamma_ebic <-gammlist[best_ebic_index[1, 1],best_ebic_index[1, 2],]
  y_pre_ebic<-X %*% gamma_ebic
  
  best_alpha_bic <- alphas[best_bic_index[1, 1]]
  best_lambda_bic <- lambdas[best_bic_index[1, 2]]
  gamma_bic <-gammlist[best_bic_index[1, 1],best_bic_index[1, 2],]
  y_pre_bic<-X %*% gamma_bic
  
  best_alpha_aic <- alphas[best_aic_index[1, 1]]
  best_lambda_aic <- lambdas[best_aic_index[1, 2]]
  gamma_aic <-gammlist[best_aic_index[1, 1],best_aic_index[1, 2],]
  y_pre_aic<-X %*% gamma_aic
  
  best_alpha_mape <- alphas[best_mape_index[1, 1]]
  best_lambda_mape <- lambdas[best_mape_index[1, 2]]
  gamma_mape <-gammlist[best_mape_index[1, 1],best_mape_index[1, 2],]
  y_pre_mape<-X %*% gamma_mape
  
  best_alpha_mae <- alphas[best_mae_index[1, 1]]
  best_lambda_mae <- lambdas[best_mae_index[1, 2]]
  gamma_mae <-gammlist[best_mae_index[1, 1],best_mae_index[1, 2],]
  y_pre_mae<-X %*% gamma_mae
  
  best_alpha_rmse <- alphas[best_rmse_index[1, 1]]
  best_lambda_rmse <- lambdas[best_rmse_index[1, 2]]
  gamma_rmse <-gammlist[best_rmse_index[1, 1],best_rmse_index[1, 2],]
  y_pre_rmse<-X %*% gamma_rmse
  
  ebic<-perform_metrics[best_ebic_index[1, 1],best_ebic_index[1, 2],1]
  bic<-perform_metrics[best_bic_index[1, 1],best_bic_index[1, 2],2]
  aic<-perform_metrics[best_aic_index[1, 1],best_aic_index[1, 2],3]
  mape<-perform_metrics[best_mape_index[1, 1],best_mape_index[1, 2],4]
  mae<-perform_metrics[best_mae_index[1, 1],best_mae_index[1, 2],5]
  rmse<-perform_metrics[best_rmse_index[1, 1],best_rmse_index[1, 2],6]
  
  result_ebic<-list(ebic,best_alpha_ebic,best_lambda_ebic,gamma_ebic,y_pre_ebic)
  names(result_ebic)<-c("ebic","alpha","lambda","gamma","y_pre")
  
  result_bic<-list(bic,best_alpha_bic,best_lambda_bic,gamma_bic,y_pre_bic)
  names(result_bic)<-c("bic","alpha","lambda","gamma","y_pre")
  
  result_aic<-list(aic,best_alpha_aic,best_lambda_aic,gamma_aic,y_pre_aic)
  names(result_aic)<-c("aic","alpha","lambda","gamma","y_pre")
  
  result_mape<-list(mape,best_alpha_mape,best_lambda_mape,gamma_mape,y_pre_mape)
  names(result_mape)<-c("mape","alpha","lambda","gamma","y_pre")
  
  result_mae<-list(mae,best_alpha_mae,best_lambda_mae,gamma_mae,y_pre_mae)
  names(result_mae)<-c("mae","alpha","lambda","gamma","y_pre")
  
  result_rmse<-list(rmse,best_alpha_rmse,best_lambda_rmse,gamma_rmse,y_pre_rmse)
  names(result_rmse)<-c("rmse","alpha","lambda","gamma","y_pre")
  
  cvresult<-list(result_ebic,result_bic,result_aic,result_mape,result_mae,result_rmse,M,GrandY,gammlist,y_pre_list,perform_metrics)
  names(cvresult)<-c("result_ebic","result_bic","result_aic","result_mape","result_mae","result_rmse","M","Grandzz","gammlist","y_pre_list","perform_metrics")
  return(cvresult)
}
clusterExport(cl,c())
clusterEvalQ(cl, library("fda"))
clusterEvalQ(cl, library(Metrics))
clusterExport(cl, list("in.mat","in.derivmat", "XMatstar", "T", "nbasis","n","p","GrandY","GrandB","bbasisMat","group"), envir=environment())

cvresultfinal<- parLapply(cl,phiseq,cvphi)

stopCluster(cl)

cves_ebic<-c()
cves_bic<-c()
cves_aic<-c()
cves_mape<-c()
cves_mae<-c()
cves_rmse<-c()
for (k in 1:length(phiseq))  
{
  cves_ebic<-c(cves_ebic,cvresultfinal[[k]]$result_ebic$ebic)#注意cvresultfinal的数据结构是在clusterExport里确定的
  cves_bic<-c(cves_bic,cvresultfinal[[k]]$result_bic$bic)
  cves_aic<-c(cves_aic,cvresultfinal[[k]]$result_aic$aic)
  cves_mape<-c(cves_mape,cvresultfinal[[k]]$result_mape$mape)
  cves_mae<-c(cves_mae,cvresultfinal[[k]]$result_mae$mae)
  cves_rmse<-c(cves_rmse,cvresultfinal[[k]]$result_rmse$rmse)
}

ind_huber<-which.min((cves_ebic))#
psi_huber<-phiseq[ind_huber]

M_huber<-cvresultfinal[[ind_huber]]$M
Grandzz<-cvresultfinal[[ind_huber]]$Grandzz

result_ebic<-cvresultfinal[[ind_huber]]$result_ebic
lambda_huber<-result_ebic$lambda
alpha_huber<-result_ebic$alpha
gamma<-result_ebic$gamma
pre_huber<-result_ebic$y_pre

indsel2<-which(gamma[-1]!=0)
Groupsel2<-group_no_intercept[indsel2]

grind2<-unique(Groupsel2)
varselected2<-grind2
namehuber<-name[varselected2+1]

result_ebic<-list(namehuber,psi_huber,alpha_huber,lambda_huber,gamma,pre_huber,M_huber)
names(result_ebic)<-c("featurename","psi","alpha","lambda","gamma","y_pre","M_huber")

###########################
ind_huber<-which.min((cves_bic))
psi_huber<-phiseq[ind_huber]
M_huber<-cvresultfinal[[ind_huber]]$M
Grandzz<-cvresultfinal[[ind_huber]]$Grandzz
result_bic<-cvresultfinal[[ind_huber]]$result_bic
lambda_huber<-result_bic$lambda
alpha_huber<-result_bic$alpha
gamma<-result_bic$gamma
pre_huber<-result_bic$y_pre
indsel2<-which(gamma[-1]!=0)
varselected2<-unique(group_no_intercept[indsel2])
namehuber<-name[varselected2+1]
result_bic<-list(namehuber,psi_huber,alpha_huber,lambda_huber,gamma,pre_huber,M_huber)
names(result_bic)<-c("featurename","psi","alpha","lambda","gamma","y_pre","M_huber")
###########################

###########################
ind_huber<-which.min((cves_aic))
psi_huber<-phiseq[ind_huber]
M_huber<-cvresultfinal[[ind_huber]]$M
Grandzz<-cvresultfinal[[ind_huber]]$Grandzz
result_aic<-cvresultfinal[[ind_huber]]$result_aic
lambda_huber<-result_aic$lambda
alpha_huber<-result_aic$alpha
gamma<-result_aic$gamma
pre_huber<-result_aic$y_pre
indsel2<-which(gamma[-1]!=0)
varselected2<-unique(group_no_intercept[indsel2])
namehuber<-name[varselected2+1]
result_aic<-list(namehuber,psi_huber,alpha_huber,lambda_huber,gamma,pre_huber,M_huber)
names(result_aic)<-c("featurename","psi","alpha","lambda","gamma","y_pre","M_huber")
###########################

###########################
ind_huber<-which.min((cves_mape))
psi_huber<-phiseq[ind_huber]
M_huber<-cvresultfinal[[ind_huber]]$M
Grandzz<-cvresultfinal[[ind_huber]]$Grandzz
result_mape<-cvresultfinal[[ind_huber]]$result_mape
lambda_huber<-result_mape$lambda
alpha_huber<-result_mape$alpha
gamma<-result_mape$gamma
pre_huber<-result_mape$y_pre
indsel2<-which(gamma[-1]!=0)
varselected2<-unique(group_no_intercept[indsel2])
namehuber<-name[varselected2+1]
result_mape<-list(namehuber,psi_huber,alpha_huber,lambda_huber,gamma,pre_huber,M_huber)
names(result_mape)<-c("featurename","psi","alpha","lambda","gamma","y_pre","M_huber")
###########################

###########################
ind_huber<-which.min((cves_mae))
psi_huber<-phiseq[ind_huber]
M_huber<-cvresultfinal[[ind_huber]]$M
Grandzz<-cvresultfinal[[ind_huber]]$Grandzz
result_mae<-cvresultfinal[[ind_huber]]$result_mae
lambda_huber<-result_mae$lambda
alpha_huber<-result_mae$alpha
gamma<-result_mae$gamma
pre_huber<-result_mae$y_pre
indsel2<-which(gamma[-1]!=0)
varselected2<-unique(group_no_intercept[indsel2])
namehuber<-name[varselected2+1]
result_mae<-list(namehuber,psi_huber,alpha_huber,lambda_huber,gamma,pre_huber,M_huber)
names(result_mae)<-c("featurename","psi","alpha","lambda","gamma","y_pre","M_huber")
###########################

###########################
ind_huber<-which.min((cves_rmse))
psi_huber<-phiseq[ind_huber]
M_huber<-cvresultfinal[[ind_huber]]$M
Grandzz<-cvresultfinal[[ind_huber]]$Grandzz
result_rmse<-cvresultfinal[[ind_huber]]$result_rmse
lambda_huber<-result_rmse$lambda
alpha_huber<-result_rmse$alpha
gamma<-result_rmse$gamma
pre_huber<-result_rmse$y_pre
indsel2<-which(gamma[-1]!=0)
varselected2<-unique(group_no_intercept[indsel2])
namehuber<-name[varselected2+1]
result_rmse<-list(namehuber,psi_huber,alpha_huber,lambda_huber,gamma,pre_huber,M_huber)
names(result_rmse)<-c("featurename","psi","alpha","lambda","gamma","y_pre","M_huber")
###########################


result_all<-list(result_ebic,result_bic,result_aic,result_mape,result_mae,result_rmse,GrandY,bbasis,bbasisint,cvresultfinal,phiseq,alphas,lambdas)
names(result_all)<-c("result_ebic","result_bic","result_aic","result_mape","result_mae","result_rmse","GrandY","bbasis","bbasisint","cvresultfinal","phiseq","alphaseq","lambdaseq")
return(result_all)
}