library(MASS)
library(mgcv)
library(refund)
library(fda)
library(parallel)
library(grpreg)
library(Matrix)
library(openxlsx)
library(tictoc)
library(ggplot2)
library(pracma)#quad
library(GA)
library(DEoptim)
library(pso)
library(Rsolnp)
library(nloptr)
library(Metrics)
library(glmnet)
library(gridExtra)
library(factoextra)
library(cluster)
library (vegan)
library(mclust)
library(data.table)
library(lpSolve)
library(Rglpk)
library(deaR)
library(iapws)
library(IAPWS95)
library(dplyr)
getwd()
setwd("D:/functional regression")

################################################################
coalfull<-read.xlsx("minute-Part234678-fulldata.xlsx")

plot(coalfull[which(coalfull$label==3),]$minutelist,coalfull[which(coalfull$label==3),]$mw,xlab='time index(minute)',ylab='load (mw)')

coaldata<-cbind(id=c(0),time=coalfull$minutelist,coalfull[,-c(1,37:43)])
coaldata$sum1701<-coaldata$m1701aml+coaldata$m1701bml#
coaldata$m1701aml<-NULL
coaldata$m1701bml<-NULL


temp<-read.xlsx('coaldata_part234678.xlsx')
temp<-temp[temp$sum1701>1,]
fdmintue<-temp$time
coaldata <- data.frame(id = temp$id, time = temp$time)

coaldata$mw<-temp$mw#负荷
coaldata$uf_te1501<-rowMeans(cbind(temp$te1501a,temp$te1501b))#强制风机空气温度

coaldata$pf_at1601<-rowMeans(cbind(temp$at1601a,temp$at1601b))#空气预热器前部烟气含氧量
coaldata$pf_at1602<-rowMeans(cbind(temp$at1602a,temp$at1602b,temp$at1602c,temp$at1602d))#空气预热器后部烟气含氧量
coaldata$pf_ft2902<-temp$ft2902#凝结水流量

coaldata$pf_pt2502<-rowMeans(cbind(temp$pt2502a,temp$pt2502b))#主蒸汽压力
coaldata$pf_rdlc04<-temp$rdlc04#主蒸汽温度
coaldata$pf_te2512<-rowMeans(cbind(temp$te2512a,temp$te2512b))#再热蒸汽温度
coaldata$pf_te3302<-rowMeans(cbind(temp$te3302a,temp$te3302b))#冷凝器入口循环水温度
coaldata$pf_te3304<-rowMeans(cbind(temp$te3304a,temp$te3304b))#冷凝器出口循环水温度

coaldata$of_ft1801<-temp$ft1801#供油流量
coaldata$of_ft2903<-rowMeans(cbind(temp$ft2903a,temp$ft2903b))#主给水
coaldata$of_m1501<-rowMeans(cbind(temp$m1501azi,temp$m1501bzi))#送风机电机电流
coaldata$of_m1602<-rowMeans(cbind(temp$m1602azi,temp$m1602bzi))#引风机电机电流
coaldata$of_m1702<-rowMeans(cbind(temp$m1702azi,temp$m1702bzi))#磨煤机主电机电流
coaldata$of_pt2807<-abs(rowMeans(cbind(temp$pt2807a,temp$pt2807b)))#冷凝器真空度
coaldata$of_x12nao<-rowMeans(cbind(temp$x12nao008,temp$x12nao012,temp$x12nao076))/1000#总风量

coaldata$sum1701<-temp$sum1701#给煤机的煤量

#########################################编号函数
find_idd<-function(m1,m2,y,the)
{
  i=3
  j=1
  start<-0
  index<-list()
  while (i <= length(y)-1) 
  {
    if ( (m1<m2) && (abs(m1-y[i])<the)) 
    {
      start<-i
      while(y[i] < m2+the)
      {
        if (y[i-2]+y[i-1]<=y[i]+y[i+1])
        {
          i<-i+1
        }
        else
        {
          i<i+1
          break
        }
      }
      if ((start<i-1)&& abs(y[i-1]-m2)<the)
      {
        index[[j]]<-c(start:(i-1))
        j<-j+1
      }
    }
    if ( (m1>m2) && (abs(m1-y[i])<the) )
    {
      start<-i
      while(y[i] > m2-the)
      {
        # if (y[i-2]+y[i-1]>=y[i]+y[i+1])
        if (!anyNA(y[c(i - 2, i - 1, i, i + 1)]) && y[i - 2] + y[i - 1] >= y[i] + y[i + 1])
        {
          i<-i+1
        }
        else
        {
          break
        }
      }
      if ((start<i-1)&& abs(y[i-1]-m2)<the)
      {
        index[[j]]<-c(start:(i-1))
        j<-j+1
      }
      else
      {
        i<-i+1
      }
    }
    else
    {
      i<-i+1
    }
  }
  return(index)
}
#########################################
k_down<-find_idd(130,100,coaldata$mw,10)
k_down[[12]]<-NULL
k_down[[16]]<-NULL
#########################################
k_up<-find_idd(100,130,coaldata$mw,10)
#########################################
coaldata_up<-coaldata
coaldata_down<-coaldata
for (i in 1:length(k_up)){
  coaldata_up$id[k_up[[i]]]<-i
}
for (i in 1:length(k_down)){
  coaldata_down$id[k_down[[i]]]<-i
}
#########################################
# k<-k_up
# coaldata<-coaldata_up
k<-k_down
coaldata<-coaldata_down
#########################################
coaldata_fullid<-coaldata
coaldata_fullid$color <- ifelse(coaldata_fullid$id == 0, "blue", "red")
coaldata_fullid$index<-c(1:length(coaldata_fullid$mw))
p_full <- ggplot(coaldata_fullid, aes(x = index, y = mw, color = color))+xlab("time") +ylab("mw")+ geom_point(size=1)+scale_color_identity()
#########################################
temp<-coaldata[coaldata$id!=0,]##########
p_single <- ggplot(temp, aes(x=c(1:length(temp$mw)),y=mw,group=id,color=id))+xlab("time") +ylab("mw")+ geom_line()
#########################################
coalseries<-temp$sum1701
coaldata$sum1701<-NULL

#########################################
coaldata<-coaldata[coaldata$id!=0,]
#########################################

xx<-c()
for ( i in 1:max(coaldata$id))
{
  start_time=(coaldata[coaldata$id==i,])$time[1]
  numb<-which(coaldata$id==i)
  xx<-c(xx,coaldata$time[numb[length(numb)]]-coaldata$time[numb[1]])
  for (j in numb)
  {
    coaldata$time[j]<-coaldata$time[j]-start_time
  }
}

#  L2 standardization
euclidean_norm <- function(x) {
  return(x / sqrt(sum(x^2)))
}
model_r2 <- function(actual, predicted) {
  r_squared <- 1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
  return(r_squared)
}

# coaldata[,-c(1,2)] <- apply(coaldata[,-c(1,2)], 2, euclidean_norm)
#MAPE
MAPE <- function(actual, forecast) {
  errors <- actual - forecast
  mape <- mean(abs(errors/actual))*100
  return(mape)
}



# standard_oredr<-which.min(xx)
standard<-min(xx)
for ( i in 1:max(coaldata$id))
{
  numb<-which(coaldata$id==i)
  for (j in numb)
  {
    coaldata$time[j]<-coaldata$time[j]/xx[i]*standard
  }
}

p_adjust <- ggplot()+xlab("time") +ylab("mw")+ geom_line(data = coaldata_down_full, aes(x=time,y=mw,group=id,color='blue'))
p_adjust_up <- ggplot() +
  xlab("time(s)") +
  ylab("loda(MW)") +
  geom_line(data = coaldata_up_full, aes(x = time, y = mw, group = id), color = 'dark gray') +
  geom_point(shape=18,color='dark gray')+
  theme_bw() +  
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = "black", size = 1),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) 

ggsave("C:/R_projects/select/result/p_adjust_up.png", plot = p_adjust_up, width = 6, height = 4, units = "in", dpi = 300)
ggsave("C:/R_projects/select/result/p_adjust_down.png", plot = p_adjust_down, width = 6, height = 4, units = "in", dpi = 300)
#########################################################################

# [1] "id"        "time"      "mw"##        "uf_te1501"##
# [5] "pf_at1601" "pf_at1602" "pf_ft2902" "pf_pt2502"
# [9] "pf_rdlc04" "pf_te2512" "pf_te3302" "pf_te3304"
# [13] "of_ft1801" "of_ft2903" "of_m1501"  "of_m1602" 
# [17] "of_m1702"  "of_pt2807" "of_x12nao"
coaldata_backup<-coaldata
#（Normalization）
# for ( i in 1:max(coaldata$id))
# {
#   numb<-which(coaldata$id==i)
#   for (j in 3:ncol(coaldata))
#   {
#     nummax<-max(coaldata[numb,j])
#     nummin<-min(coaldata[numb,j])
#     coaldata[numb,j]<-(coaldata[numb,j]-nummin)/(nummax-nummin)
#   }
# }

v1='mw'
v2='uf_te1501'
data<-coaldata[which(coaldata$id==1),]
time=data$time

#############define B-spline
splines <- create.bspline.basis(rangeval=c(0, standard), nbasis =standard+3,norder=4)
plot(splines,type="o")

#############use GCV to find lambda
gcv_set_v1<-vector()
gcv_set_v2<-vector()
num<-seq(-6,6,0.5)
for (i in num)
{
  gcv_set_v1=append(gcv_set_v1,lambda2gcv(i,time,data[,v1],fdPar(splines)))
  gcv_set_v2=append(gcv_set_v2,lambda2gcv(i,time,data[,v2],fdPar(splines)))
}
la_v1=10^(num[which.min(gcv_set_v1)])
la_v2=10^(num[which.min(gcv_set_v2)])

fdpara_v1=fdPar(splines,2,lambda = la_v1)
fdpara_v2=fdPar(splines,2,lambda = la_v2)

fd_v1 <- smooth.basis(argvals=time, data[,v1], fdpara_v1)$fd
fd_v2 <- smooth.basis(argvals=time, data[,v2], fdpara_v2)$fd


plot(time,data[,v1],xlim=c(0,20),ylim=c(90,140),col=c("blue"),ylab = 'Value')
title(paste("lamda_v1=",la_v1))
par(new=TRUE)
plot(fd_v1,xlim=c(0,20),ylim = c(90,140),col=c("red"))



plot(time,data[,v2],xlim=c(0,20),ylim = c(15,17),col=c("blue"),ylab = 'Value')
title(paste("lamda_v2=",la_v2))
par(new=TRUE)
plot(fd_v2,xlim=c(0,20),ylim = c(15,17),col=c("red"))


for (i in 2:max(coaldata$id))
{
  data<-coaldata[which(coaldata$id==i),]
  time=data$time

  gcv_set_v1<-vector()
  gcv_set_v2<-vector()
  num<-seq(-6,6,0.5)
  for (i in num)
  {
    gcv_set_v1=append(gcv_set_v1,lambda2gcv(i,time,data[,v1],fdPar(splines)))
    gcv_set_v2=append(gcv_set_v2,lambda2gcv(i,time,data[,v2],fdPar(splines)))
  }
  la_v1=10^(num[which.min(gcv_set_v1)])
  la_v2=10^(num[which.min(gcv_set_v2)])

  fdpara_v1=fdPar(splines,2,lambda = la_v1)
  fdpara_v2=fdPar(splines,2,lambda = la_v2)
  fd_v1_2 <- smooth.basis(argvals=time, data[,v1], fdpara_v1)$fd
  fd_v2_2 <- smooth.basis(argvals=time, data[,v2], fdpara_v2)$fd
  fd_v1$coefs=cbind(fd_v1$coefs,fd_v1_2$coefs)
  fd_v2$coefs=cbind(fd_v2$coefs,fd_v2_2$coefs)
}

pcaload = pca.fd(fd_v1,nharm=5)
pcatemperature=pca.fd(fd_v2,nharm=5)#



k_basis<-cbind(t(fd_v1$coefs),t(fd_v2$coefs))

kdata<-k_basis
#######################################KMA
y<-vector()
elbow<-vector()
resvect<-vector()
# ceiling(sqrt(length(k)))
x=c(2:ceiling(sqrt(length(k))))
for (K in x){
  dis <- get_dist(kdata, method='euclidean', stand = TRUE)
  res <- kmeans(dis,centers=K,nstart = 50)
  resvect=cbind(resvect,res$cluster)
  sil <- silhouette (res$cluster, dis)
  elbow=append(elbow,sum(res$withinss))
  res1<-summary(sil)
  mean_lun<-res1$avg.width
  y=append(y,mean_lun)
}
plot(x,y,type = "o",col=c("blue"))
title('blue is the silhouette coefficient for each K')
par(new=TRUE)
plot(x,elbow,type = "o",col=c("red"))
yy=which(y==max(y),arr.ind=TRUE)
print(paste('max SC is ',round(max(y),4),'best',yy+1,'clusters'))
clusternum<-resvect[,yy]
print(table(clusternum))
plot(c(1:max(coaldata$id)),clusternum,type='o',col=c("blue"))
title('blue is the cluster number for each id')



coaldata$cluster <- rep(NA, nrow(coaldata)) 
for (i in 1:max(coaldata$id)) {
  coaldata$cluster[coaldata$id == i] <- clusternum_down[i]
}


f_cluster_down_v1<-ggplot(coaldata_down_full, aes(x = time, y =mw, group = id, color = factor(cluster))) +
  geom_line() +
  labs(x = "Time (s)",
       y = "Load (MW)",
       color = "Cluster") +
  ylim(90, 140)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = "black", size = 1),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)
        ) +
  scale_color_manual(values = c("1" = "#5eb3d2", "2" = "#ff8c00"))

f_cluster_down_v2<-ggplot(coaldata_down_full, aes(x = time, y =uf_te1501, group = id, color = factor(cluster))) +
  geom_line() +
  labs(x = "Time (s)",
       y = "Temperature (℃)",
       color = "Cluster") +
  ylim(0,40)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = "black", size = 1),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) +
  scale_color_manual(values = c("1" = "#5eb3d2", "2" = "#ff8c00"))

ggsave("C:/R_projects/select/result/f_cluster_down_v1.png", plot = f_cluster_down_v1, width = 6, height = 4, units = "in", dpi = 300)
ggsave("C:/R_projects/select/result/f_cluster_down_v2.png", plot = f_cluster_down_v2, width = 6, height = 4, units = "in", dpi = 300)

#同一分类的Id画在一张图上
f_cluster_down_c2<-ggplot(coaldata[which(coaldata$cluster==2),], aes(x = time, y = mw, group = id)) +
  geom_line() +
  labs(title = "Load values for different IDs by clusters",
       x = "Time",
       y = "Load",
       color = "Cluster") +
  xlim(0, standard+2)+
  ylim(90, 140)+
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))

grid.arrange(f_cluster_down_v1,f_cluster_down_v2,f_cluster_down_c1,f_cluster_down_c2, ncol = 2)
grid.arrange(f_cluster_up_v1,f_cluster_up_v2,f_cluster_up_c1,f_cluster_up_c2, ncol = 2)


numc=2###对不同的cluster进行回归预测
##############找到对应的coal煤耗值
temp<-coaldata[coaldata$cluster==numc,]
y1<-vector()
for (i in which(clusternum_up==numc))############################这一块要注意
{
  y1<-append(y1,coalseries[which(coaldata$id==i)])
}
##############把id值重新改为从1开始
temp$cluster<-NULL
j=1
for (i in unique(temp$id)) {
  temp$id[temp$id == i] <- j
  j<-j+1
}

coaldata_up_full<-coaldata  # coaldata_up_full代表上升阶段的所有id的样本集合，最后一列是cluster
coaldata_up_cluster_2<-temp # coaldata_up_cluster_1代表上升阶段的类1的重排id后样本集合
y_up_cluster_2<-y1          # y_up_cluster_1代表上升阶段的类1的重排id后的y值

##下面的都没用了
# source('varselect_fixload.R')#调用固定负荷的版本进行对比
# standard<-max(coaldata_down_cluster_1$time)
# res_cluster_down_1<-FLCM.select(y_down_cluster_1,coaldata_down_cluster_1,nbasis1=7,nbasis2=standard+1,cvl=10)
# res_cluster_down_2<-FLCM.select(y_down_cluster_2,coaldata_down_cluster_2,nbasis1=7,nbasis2=standard+1,cvl=10)
# 
# standard<-max(coaldata_up_cluster_1$time)
# res_cluster_up_1<-FLCM.select(y_up_cluster_1,coaldata_up_cluster_1,nbasis1=7,nbasis2=standard+1,cvl=10)
# res_cluster_up_2<-FLCM.select(y_up_cluster_2,coaldata_up_cluster_2,nbasis1=7,nbasis2=standard+1,cvl=10)


com<-matrix(nrow = 8,ncol=3)
temp<-res_cluster_up_1
act<-temp$YY
com[1,]<-c(mae(act,temp$prescad),rmse(act,temp$prescad),MAPE(act,temp$prescad))
com[2,]<-c(mae(act,temp$premcp),rmse(act,temp$premcp),MAPE(act,temp$premcp))
temp<-res_cluster_up_2
act<-temp$YY
com[3,]<-c(mae(act,temp$prescad),rmse(act,temp$prescad),MAPE(act,temp$prescad))
com[4,]<-c(mae(act,temp$premcp),rmse(act,temp$premcp),MAPE(act,temp$premcp))
temp<-res_cluster_down_1
act<-temp$YY
com[5,]<-c(mae(act,temp$prescad),rmse(act,temp$prescad),MAPE(act,temp$prescad))
com[6,]<-c(mae(act,temp$premcp),rmse(act,temp$premcp),MAPE(act,temp$premcp))
temp<-res_cluster_down_2
act<-temp$YY
com[7,]<-c(mae(act,temp$prescad),rmse(act,temp$prescad),MAPE(act,temp$prescad))
com[8,]<-c(mae(act,temp$premcp),rmse(act,temp$premcp),MAPE(act,temp$premcp))

com_result_fd<-data.frame(com)
colnames(com_result_fd) <- c("MAE", "RMSE","MAPE")
rownames(com_result_fd) <- c("up_1_scad", "up_1_mcp","up_2_scad", "up_2_mcp","down_1_scad", "down_1_mcp","down_2_scad", "down_2_mcp")

# write.xlsx(com_result_fd,file = "comparison_fd_20241018.xlsx",rowNames=TRUE)

print(res_cluster_down_1$scad)
print(res_cluster_down_1$mcp)
print(res_cluster_down_2$scad)
print(res_cluster_down_2$mcp)
print(res_cluster_up_1$scad)
print(res_cluster_up_1$mcp)
print(res_cluster_up_2$scad)
print(res_cluster_up_2$mcp)


# cluster_result<-data.frame(id=c(1:154),load_up=clusternum_up,load_down=c(clusternum_down,rep(NA, 26)))
# write.xlsx(cluster_result, file="cluster_result.xlsx", rowNames=FALSE)

###########################################################作对比
timeinverse <- function(shijian) {
  if (shijian[1] == 2018) {
    interval <- (shijian[3] - 1) * 60 * 24 + shijian[4] * 60 + shijian[5]
    if (shijian[2] == 1) {
      x <- interval
    } else if (shijian[2] == 2) {
      x <- 31 * 24 * 60 + interval
    } else if (shijian[2] == 3) {
      x <- (31 + 28) * 24 * 60 + interval
    } else if (shijian[2] == 4) {
      x <- (31 + 28 + 31) * 24 * 60 + interval
    } else if (shijian[2] == 5) {
      x <- (31 + 28 + 31 + 30) * 24 * 60 + interval
    } else if (shijian[2] == 6) {
      x <- (31 + 28 + 31 + 30 + 31) * 24 * 60 + interval
    } else if (shijian[2] == 7) {
      x <- (31 + 28 + 31 + 30 + 31 + 30) * 24 * 60 + interval
    } else if (shijian[2] == 8) {
      x <- (31 + 28 + 31 + 30 + 31 + 30 + 31) * 24 * 60 + interval
    } else if (shijian[2] == 9) {
      x <- (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31) * 24 * 60 + interval
    } else if (shijian[2] == 10) {
      x <- (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30) * 24 * 60 + interval
    } else if (shijian[2] == 11) {
      x <- (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31) * 24 * 60 + interval
    } else if (shijian[2] == 12) {
      x <- (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30) * 24 * 60 + interval
    }
  }
  
  if (shijian[1] == 2019) {
    interval <- (shijian[3] - 1) * 60 * 24 + shijian[4] * 60 + shijian[5] + 365 * 24 * 60
    if (shijian[2] == 1) {
      x <- interval
    } else if (shijian[2] == 2) {
      x <- 31 * 24 * 60 + interval
    } else if (shijian[2] == 3) {
      x <- (31 + 28) * 24 * 60 + interval
    } else if (shijian[2] == 4) {
      x <- (31 + 28 + 31) * 24 * 60 + interval
    } else if (shijian[2] == 5) {
      x <- (31 + 28 + 31 + 30) * 24 * 60 + interval
    } else if (shijian[2] == 6) {
      x <- (31 + 28 + 31 + 30 + 31) * 24 * 60 + interval
    } else if (shijian[2] == 7) {
      x <- (31 + 28 + 31 + 30 + 31 + 30) * 24 * 60 + interval
    } else if (shijian[2] == 8) {
      x <- (31 + 28 + 31 + 30 + 31 + 30 + 31) * 24 * 60 + interval
    } else if (shijian[2] == 9) {
      x <- (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31) * 24 * 60 + interval
    } else if (shijian[2] == 10) {
      x <- (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30) * 24 * 60 + interval
    } else if (shijian[2] == 11) {
      x <- (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31) * 24 * 60 + interval
    } else if (shijian[2] == 12) {
      x <- (31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30) * 24 * 60 + interval
    }
  }
  
  return(x)
}

# mlresult<-read.xlsx("C:/R_projects/select/comparrison/comparison-data.xlsx")
# year <- month <- day <- hour <- minute <- minutelist <- c()
# for (j in seq_len(nrow(mlresult))) {
#   time <- unlist(strsplit(mlresult$time[j], split = "[ /:]"))
#   time <- as.integer(time)
#   
#   year <- c(year, time[1])
#   month <- c(month, time[2])
#   day <- c(day, time[3])
#   hour <- c(hour, time[4])
#   minute <- c(minute, time[5])
#   minutelist <- c(minutelist, timeinverse(time))
# }
# 
# mlresult$minutelist <- minutelist
# write.xlsx(mlresult, file="C:/R_projects/select/comparison/comparison-data-minutelist.xlsx", rowNames=FALSE)
mlresult<-read.xlsx("C:/R_projects/select/comparison/comparison-data-minutelist.xlsx")

##################################寻找对应的结果，matching_rows_up_1代表mlresult中上升阶段的第一个聚类的序号
matching_rows_up_2<-c()
for (i in which(clusternum_up==2))
{
  matching_rows_up_2<-append(matching_rows_up_2,which(mlresult$minutelist %in% fdmintue[k_up[[i]]]))
}

matching_rows_down_2<-c()
for (i in which(clusternum_down==2))
{
  matching_rows_down_2<-append(matching_rows_down_2,which(mlresult$minutelist %in% fdmintue[k_down[[i]]]))
}


ml_subset<-mlresult[matching_rows_down_2,]
act<-ml_subset$actual_coal
##
com<-matrix(nrow = 40,ncol=3)
for (i in c(4:43))
{
  pre<-ml_subset[,i]
  com[i-3,]<-c(mae(act,pre),rmse(act,pre),MAPE(act,pre))
}

com_result<-matrix(nrow = 8,ncol=3)

LRFS<-c(4,8,12,16,20)-3
for (i in c(0:3))
{
  com_result[i+1,]<-apply(com[LRFS+i,], 2, mean)
}
LR<-c(24,28,32,36,40)-3
for (i in c(0:3))
{
  com_result[i+5,]<-apply(com[LR+i,], 2, mean)
}

com_result_down_2<-data.frame(com_result)
colnames(com_result_down_2) <- c("MAE", "RMSE","MAPE")
rownames(com_result_down_2) <- c("LR-FS", "ETR-FS","NNR-FS","SVR-FS","LR", "ETR","NNR","SVR")

com_result_ml<-rbind(com_result_up_1,com_result_up_2,com_result_down_1,com_result_down_2)
# write.xlsx(com_result_ml,file = "comparison20241018.xlsx",rowNames=TRUE)


