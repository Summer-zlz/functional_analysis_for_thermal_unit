source("C:/R_projects/huberloss/result/varselect_fixload_huberloss.R")
load(file = "C:/R_projects/huberloss/result/result_huber_up1_20250213_1.RData")
load(file = "C:/R_projects/huberloss/result/result_huber_up2_20250213_1.RData")
load(file = "C:/R_projects/huberloss/result/result_huber_down1_20250213_1.RData")
load(file = "C:/R_projects/huberloss/result/result_huber_down2_20250213_1.RData")

# coaldata_up_cluster_1$of_pt2807<-abs(coaldata_up_cluster_1$of_pt2807)
# coaldata_up_cluster_2$of_pt2807<-abs(coaldata_up_cluster_2$of_pt2807)
# coaldata_down_cluster_1$of_pt2807<-abs(coaldata_down_cluster_1$of_pt2807)
# coaldata_down_cluster_2$of_pt2807<-abs(coaldata_down_cluster_2$of_pt2807)

#######################################
result_huber_down1_20250213_1<-FLCM.select.huber(y_up_cluster_1,coaldata_up_cluster_1,
                        nbasis1=7,
                        nbasis2=ceiling((standard+1)/2),
                        phiseq=10^seq(-2,10,l=15),
                        alphas=10^seq(-4,-2),
                        lambdas=seq(2,5,by=0.5))
#######################################

#######################################
# coaldata<-coaldata_up_cluster_1
# y<-y_up_cluster_1
# k<-k_up
# result<-result_huber_up1_20250213_1$result_rmse
# result$basis<-result_huber_up1_20250213_1$bbasis
# result$basis_int<-result_huber_up1_20250213_1$bbasisint

# coaldata<-coaldata_up_cluster_2
# y<-y_up_cluster_2
# k<-k_up
# result<-result_huber_up2_20250213_1$result_rmse
# result$basis<-result_huber_up2_20250213_1$bbasis
# result$basis_int<-result_huber_up2_20250213_1$bbasisint

# coaldata<-coaldata_down_cluster_1
# y<-y_down_cluster_1
# k<-k_down
# result<-result_huber_down1_20250213_1$result_rmse
# result$basis<-result_huber_down1_20250213_1$bbasis
# result$basis_int<-result_huber_down1_20250213_1$bbasisint

coaldata<-coaldata_down_cluster_2
y<-y_down_cluster_2
k<-k_down
result<-result_huber_down2_20250213_1$result_rmse
result$basis<-result_huber_down2_20250213_1$bbasis
result$basis_int<-result_huber_down2_20250213_1$bbasisint


#######################################
standard<-max(coaldata$time)
#######################################

coal_true<-y
coal_pr<-result$y_pre
test<-data.frame(coal_true=coal_true,coal_pr=coal_pr,time=coaldata$time,id=coaldata$id,index=c(1:length(coaldata$id)))

p_point<- ggplot(data =test,aes(x=index))+ylab('coalconsumption')+geom_point(aes(y=coal_pr),color='blue')+geom_point(aes(y=coal_true),color='black')+ggtitle("黑色：实际值；蓝色为预测")
p_line<- ggplot(data =test,aes(x=index))+ylab('coalconsumption')+geom_line(aes(y=coal_pr),color='blue')+geom_line(aes(y=coal_true),color='black')+ggtitle("黑色：实际值；蓝色为预测")

# p4 <- ggplot()+xlab("time") +ylab("m1702bzi")+ geom_line(data = coaldata, aes(x=time,y=m1702bzi,group=id,color=id))
# coalconsumption<-data.frame(coal=y,time=coaldata$time,id=coaldata$id)
# p5 <- ggplot()+xlab("time") +ylab("coalconsumption")+ geom_line(data = coalconsumption, aes(x=time,y=coal,group=id,color=id))
# p6 <- ggplot()+xlab("time") +ylab("cal")+ geom_line(data = calcium, aes(x=age,y=bsa,group=id,color=id))+ geom_line(data = df, aes(x=time,y=bsa,group=id,color=id))

nbasis1=7
nbasis2=ceiling((standard+1)/2)


basis<-result$basis
basis_int<-result$basis_int
inter_cons<-result$gamma[1]
names(inter_cons)<-NULL
inter<-result$gamma[2:(1+nbasis1)]
names(inter)<-NULL
gamahuber<-result$gamma[-c(1:(1+nbasis1))]
result_huber<-result$featurename

vanames<-names(coaldata)[-c(1:2)]

loca<-c()
for (i in result_huber)
{loca<-c(loca,which(vanames==i))}

b_huber<-list()
for (i in 1:length(loca))
{
  b_huber[[i]]<-c(result$M_huber%*%gamahuber[(((loca[i]-1)*nbasis2)+1):(loca[i]*nbasis2)])

}

beta<-function(t)#
{
  va_value<-c()
  for (i in 1:length(b_huber))
  {
    va_value<-c(va_value,sum(b_huber[[i]]*eval.basis(t,result$basis)))#basis
  }
  return(va_value)
}

two_int<-function(t)
{
  return(inter_cons+sum(inter*eval.basis(t,result$basis_int)))
}

x1<-c()
for (i in 1:length(coaldata$id))
{
  x1<-c(x1,two_int(coaldata$time[i])+sum(beta(coaldata$time[i])*coaldata[i,loca+2]))
}

###################################

fig<-function(k)
{
   na<-coaldata[,which(names(coaldata)==k)]
   ttemp<-data.frame(x=coaldata$time,y=na,id=coaldata$id)
   return(ggplot()+xlab("time") +ylab(k)+ geom_line(data = ttemp, aes(x=x,y=y,group=id,color=id)))
}

op_va<-c('of_ft1801','of_ft2903','of_m1501','of_m1602','of_m1702','of_pt2807','of_x12nao')
chongfu<-intersect(op_va,result_huber)
loca_op<-c()
for (i in chongfu)
{
  loca_op<-c(loca_op,which(vanames==i))
}
loca_x<-loca[!(loca %in% loca_op)]

xbasis<-create.bspline.basis(c(0,standard),nbasis2,norder = 4)
mul_X<-function(t,x_pq)
{
  p=length(loca_op)
  multX<-c()
  for (m in 1:p)
  {
    multX<-c(multX,sum(x_pq[((m-1)*(xbasis$nbasis)+1):(m*(xbasis$nbasis))]*eval.basis(t,xbasis)))
  }
  return(multX)
}

op_x<-function(id,t)
{
  Xt<-c()
  for (i in loca_op)
  {
    va<-coaldata[coaldata$id==id,i+2]
    time<-coaldata[coaldata$id==id,'time']
    
    splines <- create.bspline.basis(rangeval=c(0, standard), nbasis =standard+3,norder=4)
    gcv_set<-vector()
    num<-seq(-6,6,0.5)
    for (i in num)
    {
      gcv_set=append(gcv_set,lambda2gcv(i,time,va,fdPar(splines)))
    }
    la=10^(num[which.min(gcv_set)])
    fdpara=fdPar(splines,2,lambda = la)
    
    fd <- smooth.basis(argvals=time,va, fdpara)
    Xt<-c(Xt,predict(fd,t,deriv=0))
  }
  return(Xt)
}

fdpara<-fdPar(xbasis,2,lambda = 0.01)
ge_x<-function(id,t)
{
  Xt<-c()
  for (i in loca_x)
  {
    va<-coaldata[coaldata$id==id,i+2]
    fd <- smooth.basis(argvals=coaldata[coaldata$id==id,'time'],va, fdpara)
    Xt<-c(Xt,predict(fd,t,deriv=0))
  }
  return(Xt)
}



optima<-function(id,t,x_pq)
{
  time<-coaldata[coaldata$id==id,'time']
  mul<-c(mul_X(t,x_pq),as.numeric(coaldata[coaldata$id==id,(loca_x+2)][match(t,time),]))
  mul_ord<-c(loca_op,loca_x)
  sort_mul<-mul[order(mul_ord)]
  return(two_int(t)+sum(beta(t)*sort_mul))
}

optima_t<-function(id,t,x_pq)
{

  mul<-c(mul_X(t,x_pq),ge_x(id,t))
  mul_ord<-c(loca_op,loca_x)
  sort_mul<-mul[order(mul_ord)]
  return(two_int(t)+sum(beta(t)*sort_mul))
}

original<-function(id,t)
{

  mul<-c(op_x(id,t),ge_x(id,t))
  mul_ord<-c(loca_op,loca_x)
  sort_mul<-mul[order(mul_ord)]
  return(two_int(t)+sum(beta(t)*sort_mul))
}


fig_optim<-function(k,id)
{
  x_pq<-mtpoint[[id]]
  h=0.5
  x<-c()
  for (i in seq(0,standard,h)){
    x<-rbind(x,mul_X(i,x_pq))
  }
  x1<-c()
  for (i in seq(0,standard,h)){
    x1<-rbind(x1,op_x(id,i))
  }
  inde<-k
  ttemp<-data.frame(x=seq(0,standard,h),y1=x[,inde],y2=x1[,inde])
  return(ggplot(data = ttemp, aes(x=x))+xlab("Time (s)") +ylab("Value")
         +geom_line(aes(y=y1,color='optimal'))
         +geom_line(aes(y=y2,color='original'))
         # +ylim(0,100)
         +scale_color_manual(values = c("optimal" = "red", "original" = "blue"))

         +labs(title = paste("id is",id,", feature is",chongfu[k]))
         +theme_bw()
         +theme(panel.grid.major = element_blank(),  
                 panel.grid.minor = element_blank(), 
                 panel.border = element_rect(color = "black", size = 1),
                 axis.text = element_text(size = 14),
                 axis.title = element_text(size = 14)))
}

vanames[1]<-c("uf_load")

x3<-c()
for (t in seq(0,standard,0.1))
{
  x3<-rbind(x3,beta(t))
}
x3<-data.frame(x=seq(0,standard,0.1),x3)
figbeta<-function(num)#num=1,2,...,length(loca)
{
  x_column_name <- colnames(x3)[1]
  y_column_name <- colnames(x3)[num+1]
  return(ggplot(x3,aes_string(x = x_column_name, y = y_column_name))+ 
           geom_line(color = '#ff8c00')+
           # ylim(-0.75, 1.75)+
           labs(x="Time (mintute)",y="Value")+
           theme_bw() +  
           labs(title = paste("Feature:",vanames[loca[num]]))+
           theme(panel.grid.major = element_blank(),  
                 panel.grid.minor = element_blank(), 
                 panel.border = element_rect(color = "black", size = 1),
                 axis.text = element_text(size = 14),
                 axis.title = element_text(size = 14),
                 plot.title = element_text(hjust = 0.5))
  )
}


# windows()
# par(mfrow = c(2, 2))
# figbeta(1)
# figbeta(2)
# figbeta(3)
# figbeta(4)
# dev.off()



num_x<-7+1
num_y<-1

dea_res <- matrix(0, nrow = length(y), ncol = 1+num_x+num_y)
for (id in 1:max(coaldata$id))
{

  dea_data<-coaldata[which(coaldata$id==id),]
  input <- cbind(y[which(coaldata$id==id)],dea_data[,13:19])
  output <- dea_data$mw
  names(input)<-c("coal_quantity",names(input)[-1])
  
  num_dmu <- nrow(dea_data)  

  z <- matrix(0, nrow = num_dmu, ncol = 1+num_x+num_y)

  for(i in 1:num_dmu) {
    obj <- c(rep(1/num_x,num_x),rep(1/num_y,num_y),rep(0,num_dmu))  
    dir<-c(rep("<=",num_x),rep(">=",num_y),rep(">=",num_dmu+num_x+num_y))
    A <- matrix(0, nrow = num_x+num_y+num_dmu+num_x+num_y, ncol = num_x+num_y+num_dmu)

    A[1:(num_x+num_y),1:(num_x+num_y)]<-diag(c(input[i,],-output[i]))
    A[1:(num_x+num_y),(num_x+num_y+1):(num_x+num_y+num_dmu)] <- t(cbind(input,output))
    
    A[(num_x+num_y+1):(num_x+num_y+num_dmu),(num_x+num_y+1):(num_x+num_y+num_dmu)]<-diag(num_dmu)
    A[(num_x+num_y+num_dmu+1):(num_x+num_y+num_dmu+num_x+num_y),1:(num_x+num_y)]<-diag(num_x+num_y)
    rhs <- c(input[i,], output[i],rep(0,num_dmu+num_x+num_y))  # 右端值

    result_lp <- lp("max", obj, A, dir, rhs)
    z[i,] <- c(result_lp$solution[1:(num_x+num_y)],result_lp$objval)

  }
  dea_res[which(coaldata$id==id),]<-z
}
dea_res_single<-data.frame(id=coaldata$id,dea_res)
names(dea_res_single)<-c("id",names(input),paste0("inef_out", c(1:num_y)),"ineff_total")
# View(dea_res_single)

###################################################################################################
num_x<-7+1
num_intermed<-2
num_y<-1

dea_res_stage1 <- matrix(0, nrow = length(y), ncol = 1+num_x+num_intermed)
dea_res_stage2 <- matrix(0, nrow = length(y), ncol = 1+num_intermed+num_y)
dea_res_fullstage <- matrix(0, nrow = length(y), ncol = 1+num_x+num_y)
for (id in 1:max(coaldata$id))
{

  dea_data<-coaldata[which(coaldata$id==id),]
  input <- cbind(y[which(coaldata$id==id)],dea_data[,13:19])
  intermed <- cbind(dea_data$pf_pt2502,dea_data$pf_rdlc04)
  output <-dea_data$mw
  names(input)[1]<-c("coal_quantity")
  
  num_dmu <- nrow(dea_data) 

  z <- matrix(0, nrow = num_dmu, ncol = 1+num_x+num_intermed)

  for(i in 1:num_dmu) {
    obj <- c(rep(1/num_x,num_x),rep(1/num_intermed,num_intermed),rep(0,num_dmu)) 
    dir<-c(rep("<=",num_x),rep(">=",num_intermed),rep(">=",num_dmu+num_x+num_intermed))
    A <- matrix(0, nrow = num_x+num_intermed+num_dmu+num_x+num_intermed, ncol = num_x+num_intermed+num_dmu)
    
    A[1:(num_x+num_intermed),1:(num_x+num_intermed)]<-diag(c(input[i,],-intermed[i,]))
    A[1:(num_x+num_intermed),(num_x+num_intermed+1):(num_x+num_intermed+num_dmu)] <- t(cbind(input,intermed))
    A[(num_x+num_intermed+1):(num_x+num_intermed+num_dmu),(num_x+num_intermed+1):(num_x+num_intermed+num_dmu)]<-diag(num_dmu)
    
    A[(num_x+num_intermed+num_dmu+1):(num_x+num_intermed+num_dmu+num_x+num_intermed),1:(num_x+num_intermed)]<-diag(num_x+num_intermed)
    rhs <- c(input[i,], intermed[i,],rep(0,num_dmu+num_x+num_intermed))  # 右端值

    result_lp <- lp("max", obj, A, dir, rhs)
    z[i,] <- c(result_lp$solution[1:(num_x+num_intermed)],result_lp$objval)

  }
  dea_res_stage1[which(coaldata$id==id),]<-z
  

  z <- matrix(0, nrow = num_dmu, ncol = 1+num_intermed+num_y)

  for(i in 1:num_dmu) {
    obj <- c(rep(1/num_intermed,num_intermed),rep(1/num_y,num_y),rep(0,num_dmu))  # 目标函数系数
    dir<-c(rep("<=",num_intermed),rep(">=",num_y),rep(">=",num_dmu+num_intermed+num_y))
    A <- matrix(0, nrow = num_intermed+num_y+num_dmu+num_intermed+num_y, ncol = num_intermed+num_y+num_dmu)
    #系数矩阵A，先是X对应的无效率值，然后是Y对应的无效率值，然后是权重
    A[1:(num_intermed+num_y),1:(num_intermed+num_y)]<-diag(c(intermed[i,],-output[i]))
    A[1:(num_intermed+num_y),(num_intermed+num_y+1):(num_intermed+num_y+num_dmu)] <- t(cbind(intermed*(1+dea_res_stage1[which(coaldata$id==id),(num_x+1):(num_x+num_intermed)]),output))
    
    A[(num_intermed+num_y+1):(num_intermed+num_y+num_dmu),(num_intermed+num_y+1):(num_intermed+num_y+num_dmu)]<-diag(num_dmu)
    A[(num_intermed+num_y+num_dmu+1):(num_intermed+num_y+num_dmu+num_intermed+num_y),1:(num_intermed+num_y)]<-diag(num_intermed+num_y)
    rhs <- c(intermed[i,], output[i],rep(0,num_dmu+num_intermed+num_y))  # 右端值
    # 求解线性规划
    result_lp <- lp("max", obj, A, dir, rhs)
    z[i,] <- c(result_lp$solution[1:(num_intermed+num_y)],result_lp$objval)
    #第一列开始是投入的无效率值，倒数第二列产出的无效率值，最后一列总无效率值
  }
  dea_res_stage2[which(coaldata$id==id),]<-z
  
  dea_res_stage1<-round(dea_res_stage1,6)
  dea_res_stage2<-round(dea_res_stage2,6)
  # 初始化总体结果矩阵
  z <- matrix(0, nrow = num_dmu, ncol = 1+num_x+num_y)
  # 循环求解线性规划
  for(i in 1:num_dmu) {
    obj <- c(rep(1/num_x,num_x),rep(1/num_y,num_y),rep(0,num_dmu))  # 目标函数系数
    dir<-c(rep("<=",num_x),rep(">=",num_y),rep(">=",num_dmu+num_x+num_y))
    A <- matrix(0, nrow = num_x+num_y+num_dmu+num_x+num_y, ncol = num_x+num_y+num_dmu)
    #系数矩阵A，先是X对应的无效率值，然后是Y对应的无效率值，然后是权重
    A[1:(num_x+num_y),1:(num_x+num_y)]<-diag(c(input[i,],-output[i]))
    A[1:(num_x+num_y),(num_x+num_y+1):(num_x+num_y+num_dmu)] <- t(cbind(input*(1-dea_res_stage1[which(coaldata$id==id),1:num_x]),output*(1+dea_res_stage2[which(coaldata$id==id),(num_intermed+1):(num_intermed+num_y)])))
    
    A[(num_x+num_y+1):(num_x+num_y+num_dmu),(num_x+num_y+1):(num_x+num_y+num_dmu)]<-diag(num_dmu)
    A[(num_x+num_y+num_dmu+1):(num_x+num_y+num_dmu+num_x+num_y),1:(num_x+num_y)]<-diag(num_x+num_y)
    rhs <- c(input[i,], output[i],rep(0,num_dmu+num_x+num_y))  # 右端值
    # 求解线性规划
    result_lp <- lp("max", obj, A, dir, rhs)
    z[i,] <- c(result_lp$solution[1:(num_x+num_y)],result_lp$objval)
    #第一列开始是投入的无效率值，倒数第二列产出的无效率值，最后一列总无效率值
  }
  dea_res_fullstage[which(coaldata$id==id),]<-z
  
}

dea_res_stage<-data.frame(id=coaldata$id,dea_res_fullstage)
names(dea_res_stage)<-c("id",names(input),paste0("inef_out", c(1:num_y)),"ineff_total")
# View(dea_res_stage)

# id=85
# dea_data<-coaldata[which(coaldata$id==id),]
# input <- cbind(y[which(coaldata$id==id)],dea_data[,13:19])
# intermed <- cbind(dea_data$pf_pt2502,dea_data$pf_rdlc04)
# output <-dea_data$mw
# names(input)[1]<-c("coal_quantity")
# bijiao<-t(cbind(input*(1-dea_res_stage1[which(coaldata$id==id),1:num_x]),output*(1+dea_res_stage2[which(coaldata$id==id),(num_intermed+1):(num_intermed+num_y)])))-t(cbind(input,output))
# View(bijiao)

###########################################################获取约束范围
va_span_limit_dea<-list()
input <- cbind(y,coaldata[,(loca_op+2)])#注意这里只是提取了被选出来的那几个变量作为Input
# start_time <- Sys.time()
for (i in 1:max(coaldata$id))
{
  # print(paste("当前id=",i))
  va_span_limit_dea[[i]]<-array(0,dim=c(length(coaldata[coaldata$id==i,'time']),length(loca_op),2))
  for (j in 1:length(coaldata[coaldata$id==i,'time']))
  {
    for ( m in 1:length(loca_op))
    {
      value_ijm<-input[which(coaldata$id==i),2:(length(loca_op)+1)][j,m]
      ineff_ijm<-dea_res_stage[which(coaldata$id==i),3:(length(loca_op)+2)][j,m]
      va_span_limit_dea[[i]][j,m,]<-c(value_ijm*(1-ineff_ijm),value_ijm)
    }
  }
  end_time <- Sys.time()
  # total_time<-difftime(end_time, start_time, units = "mins")
  # print(paste("当前用时",total_time,"分钟"))
}
#########################################################################目标函数
y_point_1<-function(x_pq){
  sumy<-c()
  yi<-result$y_pre[which(coaldata$id == id)]#####改用预测值
  time<-coaldata[coaldata$id==id,'time']
  for (j in time)
  {
    sumy<-c(sumy,optima(id,j,x_pq))##注意这里其他变量用的是实际值
  }
  return(sum(sumy-yi)/sum(yi))######################注意这个目标函数
}
#########################################################################线性规划
calculate_ideal_dryness <- function(P_ms, T_ms, P_cond) {
  what<-c("s")
  #单位P_ms(MPa), T_ms(摄氏度), P_cond(MPa)=(101.325+原始值)/1000
  # 主蒸汽入口状态
  T_ms_K <- T_ms + 273.15
  P_cond=(101.325+P_cond)/1000
  s_ms <- iapws95(what, t = c(T_ms_K), p = c(P_ms)) # 主蒸汽熵 (kJ/kg·K)
  # 饱和液体和饱和蒸汽的熵
  s_f <- iapws95_sat(what, p = c(P_cond))[1]
  s_g <- iapws95_sat(what, p = c(P_cond))[2]
  # 理想排汽干度
  x_s <- (s_ms - s_f) / (s_g - s_f)
  return(x_s)
}
find_Pcond_for_threshold <- function(P_ms, T_ms, threshold,
                                     Pcond_min_raw = 0,    # 原始 P_cond 下限 (例如 0 kPa)
                                     Pcond_max_raw = 100, # 原始 P_cond 上限 (例如 50,000 kPa)
                                     nscan = 100,
                                     tol = 1e-6) 
{
  # 定义封装函数：根据原始 P_cond 调用你的 calculate_ideal_dryness
  x_from_Praw <- function(P_raw) {
    val <- tryCatch({
      calculate_ideal_dryness(P_ms = P_ms, T_ms = T_ms, P_cond = P_raw)
    }, error = function(e) NA_real_)
    return(as.numeric(val))
  }
  
  # 先粗扫：把原始 P_cond 映射成绝对压力范围大约 0.0006 ~ 22 MPa
  P_grid_raw <- seq(Pcond_min_raw, Pcond_max_raw, length.out = nscan)
  x_grid <- vapply(P_grid_raw, x_from_Praw, numeric(1))
  
  # 检查是否有有效数据
  if (all(is.na(x_grid))) {
    warning("All x_s are NA. Check input range or IAPWS function.")
    return(NA_real_)
  }
  
  # 找到第一个 x_s > threshold 的位置
  idx <- which(x_grid > threshold)
  if (length(idx) == 0) {
    return(NA_real_)  # 没有满足阈值的点
  }
  
  first_idx <- idx[1]
  if (first_idx == 1) {
    return(P_grid_raw[1])  # 最低点就满足
  }
  
  leftP_raw  <- P_grid_raw[first_idx - 1]
  rightP_raw <- P_grid_raw[first_idx]
  
  # 定义 root 函数
  froot <- function(P_raw) {
    val <- x_from_Praw(P_raw)
    if (is.na(val)) stop("x_s returned NA.")
    return(val - threshold)
  }
  
  # 检查符号变化
  fa <- froot(leftP_raw)
  fb <- froot(rightP_raw)
  if (fa > 0) return(leftP_raw)
  if (fb <= 0) return(NA_real_)
  
  # 求根（在原始单位下）
  root_raw <- tryCatch({
    uniroot(froot, lower = leftP_raw, upper = rightP_raw, tol = tol)$root
  }, error = function(e) {
    warning("uniroot failed; returning approximate grid value.")
    return(rightP_raw)
  })
  
  return(as.numeric(root_raw))
}
find_Pcond_intervals_for_threshold <- function(P_ms, T_ms, threshold,
                                               Pcond_min_raw = 0,
                                               Pcond_max_raw = 100,
                                               nscan = 1000,
                                               tol = 1e-6,
                                               merge_gap = 0.1) 
{
  # 调用原来的干度计算函数
  x_from_Praw <- function(P_raw) {
    val <- tryCatch({
      calculate_ideal_dryness(P_ms = P_ms, T_ms = T_ms, P_cond = P_raw)
    }, error = function(e) NA_real_)
    return(as.numeric(val))
  }
  
  # 粗采样
  P_grid_raw <- seq(Pcond_min_raw, Pcond_max_raw, length.out = nscan)
  x_grid <- vapply(P_grid_raw, x_from_Praw, numeric(1))
  
  if (all(is.na(x_grid))) {
    warning("All x_s are NA. Check input range or IAPWS setup.")
    return(data.frame(start_raw = numeric(0), end_raw = numeric(0),
                      start_MPa = numeric(0), end_MPa = numeric(0)))
  }
  
  # 找出超过阈值的点
  mask <- (x_grid > threshold) & !is.na(x_grid)
  if (all(!mask)) {
    message("No P_cond satisfies x_s > threshold.")
    return(data.frame(start_raw = numeric(0), end_raw = numeric(0),
                      start_MPa = numeric(0), end_MPa = numeric(0)))
  }
  
  idx <- which(mask)
  # 找出连续段
  breaks <- c(0, which(diff(idx) > 1), length(idx))
  intervals_idx <- mapply(function(i1, i2) idx[(i1 + 1):i2],
                          head(breaks, -1), tail(breaks, -1),
                          SIMPLIFY = FALSE)
  
  froot <- function(P_raw) {
    val <- x_from_Praw(P_raw)
    if (is.na(val)) stop("x_s NA")
    return(val - threshold)
  }
  
  # 逐段求边界
  raw_intervals <- lapply(intervals_idx, function(id_range) {
    left_idx <- min(id_range)
    right_idx <- max(id_range)
    left_before <- max(1, left_idx - 1)
    right_after <- min(length(P_grid_raw), right_idx + 1)
    
    P_left <- P_grid_raw[left_before]
    P_start <- P_grid_raw[left_idx]
    P_end   <- P_grid_raw[right_idx]
    P_right <- P_grid_raw[right_after]
    
    left_bound <- tryCatch({
      if (left_before < left_idx)
        uniroot(froot, lower = P_left, upper = P_start, tol = tol)$root
      else P_start
    }, error = function(e) NA_real_)
    
    right_bound <- tryCatch({
      if (right_after > right_idx)
        uniroot(froot, lower = P_end, upper = P_right, tol = tol)$root
      else P_end
    }, error = function(e) NA_real_)
    
    data.frame(start_raw = left_bound, end_raw = right_bound)
  })
  
  df <- do.call(rbind, raw_intervals)
  
  #  去除 NA 行
  df <- df[complete.cases(df), ]
  if (nrow(df) == 0) {
    warning("All intervals invalid (NA boundaries).")
    return(data.frame(start_raw = numeric(0), end_raw = numeric(0),
                      start_MPa = numeric(0), end_MPa = numeric(0)))
  }
  
  #  按起点排序
  df <- df[order(df$start_raw), ]
  
  # 合并相邻或重叠区间
  merged <- list()
  cur <- df[1, , drop = FALSE]
  for (i in 2:nrow(df)) {
    next_i <- df[i, , drop = FALSE]
    # 确保都不是NA
    if (any(is.na(next_i$start_raw), is.na(next_i$end_raw))) next
    if ((next_i$start_raw - cur$end_raw) <= merge_gap) {
      # 合并
      cur$end_raw <- max(cur$end_raw, next_i$end_raw, na.rm = TRUE)
    } else {
      merged <- append(merged, list(cur))
      cur <- next_i
    }
  }
  merged <- append(merged, list(cur))
  df_merged <- do.call(rbind, merged)
  
  #  单位换算
  df_merged$start_MPa <- (101.325 + df_merged$start_raw) / 1000
  df_merged$end_MPa   <- (101.325 + df_merged$end_raw) / 1000
  
  rownames(df_merged) <- NULL
  return(df_merged)
}

# Ramp<-c()
# for (p in 1:length(loca_op))
# {
#   max_ramp=0
#   for (id in 1: max(coaldata$id))
#   {
#     temp<-coaldata[which(coaldata$id==id),]
#     for (t in 2:length(temp$time))
#     {
#       max_ramp<-max(max_ramp,abs(temp[t,loca_op[p]+2]-temp[t-1,loca_op[p]+2])/(temp[t,2]-temp[t-1,2]))
#     }
#   }
#   Ramp<-c(Ramp,max_ramp)
# }
Ramp <- c()
for (p in 1:length(loca_op)) {
  all_ramps <- c()  # 用来保存所有 ramp 值
  
  for (id in 1:max(coaldata$id)) {
    temp <- coaldata[which(coaldata$id == id), ]
    for (t in 2:length(temp$time)) {
      ramp_val <- abs(temp[t, loca_op[p] + 2] - temp[t - 1, loca_op[p] + 2]) /
        (temp[t, 2] - temp[t - 1, 2])
      all_ramps <- c(all_ramps, ramp_val)
    }
  }
  
  # 取 99% 分位数
  q99 <- quantile(all_ramps, probs = 0.95, na.rm = TRUE)
  Ramp <- c(Ramp, q99)
}

threshold_xs<-1000#这个只是个初始值
for (id in 1:length(coaldata$id))
{
  threshold_xs=min(threshold_xs,calculate_ideal_dryness(P_ms=coaldata[id,'pf_pt2502'], T_ms=coaldata[id,'pf_rdlc04'],P_cond = coaldata[id,'of_pt2807']) )
}
for (id in 1:max(coaldata$id))
{
  print(paste("id=",id))
  time<-coaldata[coaldata$id==id,'time']
  va_span_limit_dea_and_xs<-va_span_limit_dea
  pos <- which(chongfu == "of_pt2807")
  if (length(pos) > 0) 
  {
    temp<-coaldata[which(coaldata$id==id),]
    for (t in 1:length(time))
    {
      p_cond_for_x_s<-find_Pcond_intervals_for_threshold(P_ms=temp[t,'pf_pt2502'], T_ms=temp[t,'pf_rdlc04'],threshold=threshold_xs)[1,1]
      va_span_limit_dea_and_xs[[id]][t,pos,1]<-max(p_cond_for_x_s,va_span_limit_dea[[id]][t,pos,1])
    }
  } 
}
# set.seed(42)
# delta<-abs(rnorm(length(loca_op)*(xbasis$nbasis), mean = 0, sd = 0.01))

delta_ratio<-c(0.00,0.10,0.20)
gamma_robust<-c(1.0,0.8,0.6)

delta_gamma_matrix<-matrix(0,nrow=length(delta_ratio),ncol=length(gamma_robust))
delta_gamma_matrix_true<-delta_gamma_matrix
for (d in 1:length(delta_ratio))
{
  for (g in 1:length(gamma_robust))
  {
    dea_lp_modify<-function(id)
    {
      delta<-c()
      for (p in 1:length(loca_op))
      {
        delta<-c(delta,abs(b_huber[[p]]*delta_ratio[d]))
      }
      
      time<-coaldata[coaldata$id==id,'time']
      
      fobj<-c()
      for (p in match(loca_op, loca))
      {
        tttt<-c(rep(0,xbasis$nbasis))
        for (t in time)
        {
          tttt<-tttt+beta(t)[p]*eval.basis(t,xbasis)
        }
        fobj<-c(fobj,tttt)
      }
      fobj<-c(fobj,rep(0,length(time)*length(loca_op)),delta)
      
      ncol =length(loca_op)*(xbasis$nbasis)+length(time)*length(loca_op)+length(loca_op)*(xbasis$nbasis)
      nrow = length(time)*length(loca_op)
      fcons_ub<-matrix(0, nrow = nrow, ncol =ncol)
      for (p in 1:length(loca_op))
      {
        for (t in 1:length(time))
        {
          fcons_ub[(p-1)*length(time)+t,(1+(p-1)*(xbasis$nbasis)) : (p*(xbasis$nbasis))]<-c(eval.basis(time[t],xbasis))
        }
      }
      
      #顺序：变量1的各个时间点，变量2的各个时间点
      frhs_ub<-c()
      for (p in c(1:length(loca_op)))
      {
        frhs_ub<-c(frhs_ub,va_span_limit_dea_and_xs[[id]][,p,2])
      }
      
      frhs_lb<-c()
      for (p in c(1:length(loca_op)))
      {
        frhs_lb<-c(frhs_lb,va_span_limit_dea_and_xs[[id]][,p,1])
      }
      
      #对于mu来说，p=1，然后time个时间，之后p=2，循环
      nrow = length(loca_op)
      
      fcon_mu_gamma<-matrix(0, nrow = nrow, ncol =ncol)
      frhs_mu_gamma<-rep(0,length(loca_op))
      for (p in 1:length(loca_op))
      {
        fcon_mu_gamma[p,(length(loca_op)*(xbasis$nbasis)+(p-1)*length(time)+1):(length(loca_op)*(xbasis$nbasis)+(p)*length(time))]<-1
        frhs_mu_gamma[p]<-floor(gamma_robust[g]*length(time))
      }
      
      nrow = length(time)*length(loca_op)
      
      fcon_mu_lb<-matrix(0, nrow = nrow, ncol =ncol)
      for (p in 1:length(loca_op))
      {
        for (t in 1:length(time))
        {
          fcon_mu_lb[(p-1)*length(time)+t,(1+(p-1)*(xbasis$nbasis)) : (p*(xbasis$nbasis))]<-c(eval.basis(time[t],xbasis))
        }
        fcon_mu_lb[((p-1)*length(time)+1):(p*length(time)),(length(loca_op)*(xbasis$nbasis)+(p-1)*length(time)+1):(length(loca_op)*(xbasis$nbasis)+(p)*length(time))]<-diag(frhs_ub[((p-1)*length(time)+1):(p*length(time))]-frhs_lb[((p-1)*length(time)+1):(p*length(time))])
      }
      frhs_mu_lb<-frhs_ub
      
      ceta<-matrix(0,nrow=xbasis$nbasis,ncol=1)
      Bp<-matrix(0,nrow=xbasis$nbasis,ncol=xbasis$nbasis)
      for (t in 1:length(time))
      {
        ceta[,1]<-c(eval.basis(time[t],xbasis))
        Bp<-Bp+ceta%*%t(ceta)
      }
      nrow = length(loca_op)*(xbasis$nbasis)
      fcon_B_S_lb<-matrix(0, nrow = nrow, ncol =ncol)
      for (p in 1:length(loca_op))
      {
        fcon_B_S_lb[((p-1)*(xbasis$nbasis)+1):((p)*(xbasis$nbasis)),((p-1)*(xbasis$nbasis)+1):((p)*(xbasis$nbasis))]<-t(Bp)
        fcon_B_S_lb[((p-1)*(xbasis$nbasis)+1):((p)*(xbasis$nbasis)),(length(loca_op)*(xbasis$nbasis)+length(time)*length(loca_op)+(p-1)*(xbasis$nbasis)+1):(length(loca_op)*(xbasis$nbasis)+length(time)*length(loca_op)+(p)*(xbasis$nbasis))]<-diag(1,xbasis$nbasis)
      }
      fcon_B_S_ub<-matrix(0, nrow = nrow, ncol =ncol)
      for (p in 1:length(loca_op))
      {
        fcon_B_S_ub[((p-1)*(xbasis$nbasis)+1):((p)*(xbasis$nbasis)),((p-1)*(xbasis$nbasis)+1):((p)*(xbasis$nbasis))]<-t(Bp)
        fcon_B_S_ub[((p-1)*(xbasis$nbasis)+1):((p)*(xbasis$nbasis)),(length(loca_op)*(xbasis$nbasis)+length(time)*length(loca_op)+(p-1)*(xbasis$nbasis)+1):(length(loca_op)*(xbasis$nbasis)+length(time)*length(loca_op)+(p)*(xbasis$nbasis))]<-diag(-1,xbasis$nbasis)
      }
      frhs_B_S_lb<-rep(0,length(loca_op)*(xbasis$nbasis))
      frhs_B_S_ub<-rep(0,length(loca_op)*(xbasis$nbasis))
      
      nrow = length(loca_op)*(length(time)-1)
      
      fcon_ramp<-matrix(0, nrow = nrow, ncol =ncol)
      frhs_ramp_down<-c()
      frhs_ramp_up<-c()
      for (p in 1:length(loca_op))
      {
        for (t in 2:length(time))
        {
          fcon_ramp[(p-1)*(length(time)-1)+t-1,(1+(p-1)*(xbasis$nbasis)) : (p*(xbasis$nbasis))]<-c(eval.basis(time[t],xbasis))-c(eval.basis(time[t-1],xbasis))
          frhs_ramp_down<-c(frhs_ramp_down,-Ramp[p]*(time[t]-time[t-1]))
          frhs_ramp_up<-c(frhs_ramp_up,Ramp[p]*(time[t]-time[t-1]))
        }
        # frhs_ramp_down<-c(frhs_ramp_down,rep(-Ramp[p],length(time)-1))
        # frhs_ramp_up<-c(frhs_ramp_up,rep(Ramp[p],length(time)-1))
      }
      
      fcon<-rbind(fcon_mu_gamma,fcons_ub,fcon_mu_lb,fcon_B_S_lb,fcon_B_S_ub,fcon_ramp,fcon_ramp)
      frhs<-c(frhs_mu_gamma,frhs_ub,frhs_mu_lb,frhs_B_S_lb,frhs_B_S_ub,frhs_ramp_down,frhs_ramp_up)
      
      fdir<-c(rep("<=",length(loca_op)+length(time)*length(loca_op)),rep(">=",length(time)*length(loca_op)),rep(">=",length(loca_op)*(xbasis$nbasis)),rep("<=",length(loca_op)*(xbasis$nbasis)),rep(">=",length(loca_op)*(length(time)-1)),rep("<=",length(loca_op)*(length(time)-1)))
      
      # res <- lp("min", fobj, fcon, fdir, frhs,binary.vec=c((length(loca_op)*(xbasis$nbasis)+1):(length(loca_op)*(xbasis$nbasis)+length(time)*length(loca_op))))
      types<-c(rep("C",length(loca_op)*(xbasis$nbasis)),rep("B",length(time)*length(loca_op)),rep("C",length(loca_op)*(xbasis$nbasis)))
      res <- Rglpk_solve_LP( fobj, fcon, fdir, frhs,types=types,max = FALSE,control = list(tm_limit = 60000))
      
      return(res)
    }
    vpoint<-c()
    mtpoint<-list()
    lastpoint<-list()
    idnum<-c(1:max(coaldata$id))
    for (id in idnum)
    {
      res<-dea_lp_modify(id)
      print(res$status)
      if(res$status==0){
          vpoint[id]<-y_point_1(res$solution)
          print(paste("id=",id,",reduction potential:",y_point_1(res$solution)))
      }
      else
      {
          vpoint[id]<-1
          # print("lp is fail")
      }
      lastpoint[[id]]<-res
      mtpoint[[id]]<-res$solution
    }
    delta_gamma_matrix[d,g]<-mean(vpoint[vpoint<0])
    delta_gamma_matrix_true[d,g]<-table(vpoint<0)[2]
    print(paste("delta=",delta_ratio[d],"gamma rate=",gamma_robust[g],"average reduction",delta_gamma_matrix[d,g],"achieved amount",delta_gamma_matrix_true[d,g]))
  }
}
write.xlsx(data.frame(delta_gamma_matrix),file='delta_gamma_matrix_down_2_1211.xlsx')


dea_lp_modify<-function(id)
{
  delta<-c()
  for (p in 1:length(loca_op))
  {
    delta<-c(delta,abs(b_huber[[p]]*delta_ratio))
  }
  
  time<-coaldata[coaldata$id==id,'time']
  
  va_span_limit_dea_and_xs<-va_span_limit_dea
  pos <- which(loca_op == 16)
  if (length(pos) > 0) 
  {
    temp<-coaldata[which(coaldata$id==id),]
    for (t in 1:length(time))
    {
      p_cond_for_x_s<-find_Pcond_intervals_for_threshold(P_ms=temp[t,'pf_pt2502'], T_ms=temp[t,'pf_rdlc04'],threshold=threshold_xs)[1,1]
      va_span_limit_dea_and_xs[[id]][t,pos,1]<-max(p_cond_for_x_s,va_span_limit_dea[[id]][t,pos,1])
    }
  } 
  
  fobj<-c()
  for (p in match(loca_op, loca))
  {
    tttt<-c(rep(0,xbasis$nbasis))
    for (t in time)
    {
      tttt<-tttt+beta(t)[p]*eval.basis(t,xbasis)
    }
    fobj<-c(fobj,tttt)
  }
  fobj<-c(fobj,rep(0,length(time)*length(loca_op)),delta)
  
  ncol =length(loca_op)*(xbasis$nbasis)+length(time)*length(loca_op)+length(loca_op)*(xbasis$nbasis)
  nrow = length(time)*length(loca_op)
  fcons_ub<-matrix(0, nrow = nrow, ncol =ncol)
  for (p in 1:length(loca_op))
  {
    for (t in 1:length(time))
    {
      fcons_ub[(p-1)*length(time)+t,(1+(p-1)*(xbasis$nbasis)) : (p*(xbasis$nbasis))]<-c(eval.basis(time[t],xbasis))
    }
  }
  
  #顺序：变量1的各个时间点，变量2的各个时间点
  frhs_ub<-c()
  for (p in c(1:length(loca_op)))
  {
    frhs_ub<-c(frhs_ub,va_span_limit_dea_and_xs[[id]][,p,2])
  }
  
  frhs_lb<-c()
  for (p in c(1:length(loca_op)))
  {
    frhs_lb<-c(frhs_lb,va_span_limit_dea_and_xs[[id]][,p,1])
  }
  
  #对于mu来说，p=1，然后time个时间，之后p=2，循环
  nrow = length(loca_op)

  fcon_mu_gamma<-matrix(0, nrow = nrow, ncol =ncol)
  frhs_mu_gamma<-rep(0,length(loca_op))
  for (p in 1:length(loca_op))
  {
    fcon_mu_gamma[p,(length(loca_op)*(xbasis$nbasis)+(p-1)*length(time)+1):(length(loca_op)*(xbasis$nbasis)+(p)*length(time))]<-1
    frhs_mu_gamma[p]<-floor(gamma_robust*length(time))
  }
  
  nrow = length(time)*length(loca_op)

  fcon_mu_lb<-matrix(0, nrow = nrow, ncol =ncol)
  for (p in 1:length(loca_op))
  {
    for (t in 1:length(time))
    {
     fcon_mu_lb[(p-1)*length(time)+t,(1+(p-1)*(xbasis$nbasis)) : (p*(xbasis$nbasis))]<-c(eval.basis(time[t],xbasis))
    }
    fcon_mu_lb[((p-1)*length(time)+1):(p*length(time)),(length(loca_op)*(xbasis$nbasis)+(p-1)*length(time)+1):(length(loca_op)*(xbasis$nbasis)+(p)*length(time))]<-diag(frhs_ub[((p-1)*length(time)+1):(p*length(time))]-frhs_lb[((p-1)*length(time)+1):(p*length(time))])
  }
  frhs_mu_lb<-frhs_ub

  ceta<-matrix(0,nrow=xbasis$nbasis,ncol=1)
  Bp<-matrix(0,nrow=xbasis$nbasis,ncol=xbasis$nbasis)
  for (t in 1:length(time))
  {
    ceta[,1]<-c(eval.basis(time[t],xbasis))
    Bp<-Bp+ceta%*%t(ceta)
  }
  nrow = length(loca_op)*(xbasis$nbasis)
  fcon_B_S_lb<-matrix(0, nrow = nrow, ncol =ncol)
  for (p in 1:length(loca_op))
  {
    fcon_B_S_lb[((p-1)*(xbasis$nbasis)+1):((p)*(xbasis$nbasis)),((p-1)*(xbasis$nbasis)+1):((p)*(xbasis$nbasis))]<-t(Bp)
    fcon_B_S_lb[((p-1)*(xbasis$nbasis)+1):((p)*(xbasis$nbasis)),(length(loca_op)*(xbasis$nbasis)+length(time)*length(loca_op)+(p-1)*(xbasis$nbasis)+1):(length(loca_op)*(xbasis$nbasis)+length(time)*length(loca_op)+(p)*(xbasis$nbasis))]<-diag(1,xbasis$nbasis)
  }
  fcon_B_S_ub<-matrix(0, nrow = nrow, ncol =ncol)
  for (p in 1:length(loca_op))
  {
    fcon_B_S_ub[((p-1)*(xbasis$nbasis)+1):((p)*(xbasis$nbasis)),((p-1)*(xbasis$nbasis)+1):((p)*(xbasis$nbasis))]<-t(Bp)
    fcon_B_S_ub[((p-1)*(xbasis$nbasis)+1):((p)*(xbasis$nbasis)),(length(loca_op)*(xbasis$nbasis)+length(time)*length(loca_op)+(p-1)*(xbasis$nbasis)+1):(length(loca_op)*(xbasis$nbasis)+length(time)*length(loca_op)+(p)*(xbasis$nbasis))]<-diag(-1,xbasis$nbasis)
  }
  frhs_B_S_lb<-rep(0,length(loca_op)*(xbasis$nbasis))
  frhs_B_S_ub<-rep(0,length(loca_op)*(xbasis$nbasis))
  
  nrow = length(loca_op)*(length(time)-1)

  fcon_ramp<-matrix(0, nrow = nrow, ncol =ncol)
  frhs_ramp_down<-c()
  frhs_ramp_up<-c()
  for (p in 1:length(loca_op))
  {
    for (t in 2:length(time))
    {
      fcon_ramp[(p-1)*(length(time)-1)+t-1,(1+(p-1)*(xbasis$nbasis)) : (p*(xbasis$nbasis))]<-c(eval.basis(time[t],xbasis))-c(eval.basis(time[t-1],xbasis))
      frhs_ramp_down<-c(frhs_ramp_down,-Ramp[p]*(time[t]-time[t-1]))
      frhs_ramp_up<-c(frhs_ramp_up,Ramp[p]*(time[t]-time[t-1]))
    }
    # frhs_ramp_down<-c(frhs_ramp_down,rep(-Ramp[p],length(time)-1))
    # frhs_ramp_up<-c(frhs_ramp_up,rep(Ramp[p],length(time)-1))
  }

  fcon<-rbind(fcon_mu_gamma,fcons_ub,fcon_mu_lb,fcon_B_S_lb,fcon_B_S_ub,fcon_ramp,fcon_ramp)
  frhs<-c(frhs_mu_gamma,frhs_ub,frhs_mu_lb,frhs_B_S_lb,frhs_B_S_ub,frhs_ramp_down,frhs_ramp_up)
  
  fdir<-c(rep("<=",length(loca_op)+length(time)*length(loca_op)),rep(">=",length(time)*length(loca_op)),rep(">=",length(loca_op)*(xbasis$nbasis)),rep("<=",length(loca_op)*(xbasis$nbasis)),rep(">=",length(loca_op)*(length(time)-1)),rep("<=",length(loca_op)*(length(time)-1)))
  
  res <- lp("min", fobj, fcon, fdir, frhs,binary.vec=c((length(loca_op)*(xbasis$nbasis)+1):(length(loca_op)*(xbasis$nbasis)+length(time)*length(loca_op))))
  
  return(res)
}



#########################################################################画优化对比的图像
# for (k in 1:length(loca_op))
# {
#   for (i in idnum[vpoint<0][4])
#   {
#     ggsave(paste0("C:/R_projects/huberloss/result/optimal_figure/optimal_up_1_id_",i,"_",vanames[loca_op][k],"_2.png"), plot = fig_optim(k,i), width = 6, height = 4, units = "in", dpi = 300)
#   }
# }

#_1是加变量上限，_2是增加基函数
#########################################################################
# va_span_limit_up_1<-va_span_limit_dea_and_xs
# lastpoint_up_1<-lastpoint
# vpoint_up_1<-vpoint
# mtpoint_up_1<-mtpoint

# va_span_limit_up_2<-va_span_limit_dea_and_xs
# lastpoint_up_2<-lastpoint
# vpoint_up_2<-vpoint
# mtpoint_up_2<-mtpoint

# va_span_limit_down_1<-va_span_limit_dea_and_xs
# lastpoint_down_1<-lastpoint
# vpoint_down_1<-vpoint
# mtpoint_down_1<-mtpoint

va_span_limit_down_2<-va_span_limit_dea_and_xs
lastpoint_down_2<-lastpoint
vpoint_down_2<-vpoint
mtpoint_down_2<-mtpoint


# #########################################################################
# z<-data.frame(id=c(1:length(vpoint_up_1)),reduction_rate=vpoint_up_1)
# write.xlsx(z,file = "C:/R_projects/huberloss/result/reduction_rate_up_1_mdea.xlsx",rowNames=FALSE)
# z<-data.frame(id=c(1:length(vpoint_up_2)),reduction_rate=vpoint_up_2)
# write.xlsx(z,file = "C:/R_projects/huberloss/result/reduction_rate_up_2_mdea.xlsx",rowNames=FALSE)
# z<-data.frame(id=c(1:length(vpoint_down_1)),reduction_rate=vpoint_down_1)
# write.xlsx(z,file = "C:/R_projects/huberloss/result/reduction_rate_down_1_mdea.xlsx",rowNames=FALSE)
# z<-data.frame(id=c(1:length(vpoint_down_2)),reduction_rate=vpoint_down_2)
# write.xlsx(z,file = "C:/R_projects/huberloss/result/reduction_rate_down_2_mdea.xlsx",rowNames=FALSE)


# mean(vpoint_up_1[vpoint_up_1<0])
# table(vpoint_up_1<0)
# mean(vpoint_up_2[vpoint_up_2<0])
# table(vpoint_up_2<0)
# mean(vpoint_down_1[vpoint_down_1<0])
# table(vpoint_down_1<0)
# mean(vpoint_down_2[vpoint_down_2<0])
# table(vpoint_down_2<0)
# 
# y_preee<-result_huber_up1_20250213_1$result_rmse$y_pre
# y_actual<-y_up_cluster_1
# mae(y_actual,y_preee)
# rmse(y_actual,y_preee)
# mape(y_actual,y_preee)
# 
# y_preee<-result_huber_up2_20250213_1$result_rmse$y_pre
# y_actual<-y_up_cluster_2
# mae(y_actual,y_preee)
# rmse(y_actual,y_preee)
# mape(y_actual,y_preee)
# 
# y_preee<-result_huber_down1_20250213_1$result_rmse$y_pre
# y_actual<-y_down_cluster_1
# mae(y_actual,y_preee)
# rmse(y_actual,y_preee)
# mape(y_actual,y_preee)
# 
# y_preee<-result_huber_down2_20250213_1$result_rmse$y_pre
# y_actual<-y_down_cluster_2
# mae(y_actual,y_preee)
# rmse(y_actual,y_preee)
# mape(y_actual,y_preee)



optima_op_data<-list()
for (i in 1:max(coaldata$id))
{
  time<-coaldata[coaldata$id==i,'time']
  optima_op_data[[i]]<-array(0,dim=c(length(time),length(loca_op)))
  for (j in 1:length(time))
  {
    for ( p in 1:length(loca_op))
    {
      optima_op_data[[i]][j,p]<-sum(mtpoint[[i]][(1+(p-1)*(xbasis$nbasis)) : (p*(xbasis$nbasis))]*c(eval.basis(time[j],xbasis)))
    }
  }
}
#下面是导出用excel画
# write.xlsx(data.frame(va_span_limit_dea[[1]]),file = "C:/R_projects/huberloss/result/up1_id1_op1_limit.xlsx",rowNames=FALSE)
# write.xlsx(data.frame(optima_op_data[[1]]),file = "C:/R_projects/huberloss/result/up1_id1_op1_optimal.xlsx",rowNames=FALSE)

last_figure<-function(id,p){
  LB<-va_span_limit_dea_and_xs[[id]][,p,1]
  UB<-va_span_limit_dea_and_xs[[id]][,p,2]
  Adjusted<-optima_op_data[[id]][,p]
  figuredata<-data.frame(x=c(1:length(LB)),LB=LB,UB=UB,Adjusted=Adjusted)
  # 将数据转换为长格式
  figuredata_long <- tidyr::pivot_longer(figuredata, cols = c(LB, UB, Adjusted), 
                                         names_to = "variable", values_to = "value")
  # 绘制图形
  p <- ggplot(data = figuredata_long, aes(x = x, y = value, color = variable, shape = variable)) +
    geom_point(size = 4) +  # 添加散点
    scale_color_manual(values = c("LB" = "blue", "UB" = "black", "Adjusted" = "red")) +  # 自定义颜色
    scale_shape_manual(values = c("LB" = 0, "UB" = 5, "Adjusted" = 8)) +  # 自定义形状
    labs(x = "Points", y = "Value", color = "Variable", shape = "Variable",title = paste("Feature:",vanames[loca_op[p]])) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),  # 隐藏主要网格线
          panel.grid.minor = element_blank(),  # 隐藏次要网格线
          panel.border = element_rect(color = "black", size = 1),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5) )+
    scale_x_continuous(labels = function(x) as.character(as.integer(x))) 
  return(p)
}
#2,6,16,https://ggplot2.tidyverse.org/reference/scale_shape.html


for (i in idnum[vpoint<0])
{
  for (j in 1:length(loca_op))
  {
    ji<-which(c(11:17)==loca_op[j])
    # ggsave(paste0("C:/R_projects/huberloss - 副本/result/optimal_figure_point/up_1/optimal_figure_up_1_id_",i,"_of_",ji,".png"), plot = last_figure(i,j), width = 6, height = 2.5, units = "in", dpi = 300)
    # ggsave(paste0("C:/R_projects/huberloss - 副本/result/optimal_figure_point/up_2/optimal_figure_up_2_id_",i,"_of_",ji,".png"), plot = last_figure(i,j), width = 6, height = 2.5, units = "in", dpi = 300)
    # ggsave(paste0("C:/R_projects/huberloss - 副本/result/optimal_figure_point/down_1/optimal_figure_down_1_id_",i,"_of_",ji,".png"), plot = last_figure(i,j), width = 6, height = 2.5, units = "in", dpi = 300)
    ggsave(paste0("C:/R_projects/huberloss - 副本/result/optimal_figure_point/down_2/optimal_figure_down_2_id_",i,"_of_",ji,".png"), plot = last_figure(i,j), width = 6, height = 2.5, units = "in", dpi = 300)
  }
}


for (i in 1:length(loca))
{
  # ggsave(paste0("C:/R_projects/huberloss/result/beta_figure/beta_up_1_",vanames[loca[i]],".png"), plot = figbeta(i), width = 6, height = 2.5, units = "in", dpi = 300)
  # ggsave(paste0("C:/R_projects/huberloss/result/beta_figure/beta_up_2_",vanames[loca[i]],".png"), plot = figbeta(i), width = 6, height = 2.5, units = "in", dpi = 300)
  # ggsave(paste0("C:/R_projects/huberloss/result/beta_figure/beta_down_1_",vanames[loca[i]],".png"), plot = figbeta(i), width = 6, height = 2.5, units = "in", dpi = 300)
  ggsave(paste0("C:/R_projects/huberloss/result/beta_figure/beta_down_2_",vanames[loca[i]],".png"), plot = figbeta(i), width = 6, height = 2.5, units = "in", dpi = 300)
}

