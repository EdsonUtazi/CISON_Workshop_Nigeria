#libraries
library(maptools); library(spdep); library(INLA)
library(mvtnorm); library(msm)
library(MASS) 

#X and E
#setwd("G:/MCMC_SOU_2/Ohio_data")
setwd("/home/edson/MCMC_SOU_2/Ohio_data")
n=88; T=21
data=read.csv("OhioRespMort.csv", header=T)
Ohio.gen = readShapePoly("tl_2010_39_county00")
hh=data.frame(Ohio.gen)
NN = hh$NAME00
NN.order = order(order(NN))

#Re-number counties to match numbering using NAME00 in shapefile
county.1=0
for (i in 1:n) county.1[i]=which(NN.order==i)

county.shp=county.shp1=data$county
for (i in 1:n) county.shp[county.shp1==i]=county.1[i]

data$county.shp=county.shp #Contains numbering used in data. SO site-specific param follow these numbering
#View(data)                 #and county names

X=matrix(0,T,n)
E=matrix(0,T,n)
for (i in 1:T)
{
  for (j in 1:n){
    X[i,j]=as.numeric(subset(data, county.shp==j & year==i, select=y))
    E[i,j]=as.numeric(subset(data, county.shp==j & year==i, select=E))
  }
}

#z
p=2
z=matrix(0,ncol=n*p,nrow=T)
a=matrix(0,n,p)
for (t in 1:T){
  for (j in 1:n) a[j,]=c(1, t/T) #change this with changing p
  z[t,]=c(t(a))}

#d
#Transform coordinates to UTM
LongLatToUTM<-function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
  return(as.data.frame(res))
}

Ohio.gen = readShapePoly("tl_2010_39_county00")
coo = coordinates(Ohio.gen)  #Check whether this is Northing/Easting. Convert if not.

#dd = LongLatToUTM(coo[,1],coo[,2],17)
#d=dist(dd, method="euclidean")  #distance in meters. divide by 1000 to get dist in km
d=dist(coo, method="euclidean")
d=as.matrix(d)

#ll
ll=matrix(0,ncol=n,nrow=T)



load("/home/edson/MCMC_SOU_3Ohio/GDDWoutcb_forbeta.rda")    #CHANGE ME
#load("/home/edson/MCMC_SOU_3Ohio/test/GDDWoutcb_forbeta.rda")

#Take average of parameters
#output[jj,]=c(1,0,taus1,alpha1,oo,tt,bb[,1:p])
#output[jj,]=c(2,theta,taus2,alpha2,oo,tt,bb[,1:p])
size = 5000                #CHANGE THIS
burn.in = (50/100)*size
rr = (burn.in + 1):size
x=X

#Predicted X, RR and prob(RR>1) - p.RR
x.pred = x.pred.out$x.pred; RR = x.pred.out$RR
x.pred = x.pred[rr,,]; RR = RR[rr,,]

x.pred.mean=RR.mean=p.RR=SMR=matrix(0,T,n)
for (t in 1:T) 
{
  uu=x.pred[,t,]; uu1=RR[,t,]
  x.pred.mean[t,]=colMeans(uu); RR.mean[t,]=colMeans(uu1)
  SMR[t,]=X[t,]/E[t,]
  for (i in 1:n){
    vv=RR[,t,i]
    p.RR[t,i] = length(which(vv>1))/(size-burn.in)
  }
}

colnames(x.pred.mean) = NN
colnames(RR.mean) = NN
colnames(p.RR) = NN
colnames(SMR) = NN


aa=data$county.shp[1:n]
x.pred.mean1 = x.pred.mean2=x.pred.mean
RR.mean1 = RR.mean2=RR.mean
p.RR1 = p.RR2 = p.RR    #probability
SMR1 = SMR2 = SMR
nname=numeric(n)
for (i in 1:n){
  x.pred.mean1[,aa[i]]=x.pred.mean2[,i]; RR.mean1[,aa[i]]=RR.mean2[,i]
  p.RR1[,aa[i]]=p.RR2[,i]
  SMR1[,aa[i]]=SMR2[,i]
  nname[aa[i]]=colnames(x.pred.mean)[i]
}
colnames(x.pred.mean1)=nname
colnames(RR.mean1)=nname
colnames(p.RR1)=nname
colnames(SMR1)=nname

#Transpose
x.pred.mean1.t = data.frame(t(x.pred.mean1)); x.pred.mean1.t$NAME00 = as.character(rownames(x.pred.mean1.t)) 
RR.mean1.t = t(RR.mean1)
p.RR1.t = t(p.RR1)
SMR1.t = t(SMR1)
SMR_ovmean = apply(SMR1.t,1,"mean")
RR_ovmean = apply(RR.mean1.t,1,"mean")
p.RR_ovmean = apply(p.RR1.t,1,"mean")

out = data.frame(x.pred.mean1.t,RR.mean1.t,p.RR1.t,SMR1.t,SMR_ovmean,RR_ovmean,p.RR_ovmean)



########################################################
#For saving and plotting
#Names
name1=name2=name3=name4=nb=0; nr=1967
for (i in 1:T){
  nb[i]=nr+i
  name1[i]=paste0("pred_",nb[i])
  name2[i]=paste0("RR_",nb[i])
  name3[i]=paste0("Prob_",nb[i])  #p.RR prob of RR>1
  name4[i]=paste0("SMR_",nb[i])
}
names=c(name1,"NAME00",name2,name3,name4,"SMRom","RRom","p.RRom")
colnames(out) = names

Ohio.gen1=merge(Ohio.gen,out, by="NAME00")

#writePolyShape(Ohio.gen1, "/home/edson/MCMC_SOU_3Ohio/OhioGDDW.shp")

library(RColorBrewer)
#SMR
pdf("/home/edson/MCMC_SOU_3Ohio/SMR.pdf") #name4[c(2,8,15,21)]
att=c(0,0.2,0.5,0.8,1.0,1.2,1.5,2.6)
spplot(Ohio.gen1, name4[c(2,8,15,21)], col.regions = brewer.pal(9, "Blues") , at=att,
       as.table=TRUE, main = "",colorkey = list(labels = list( labels = att,        
                                        at = att)))
dev.off()

#RR
pdf("/home/edson/MCMC_SOU_3Ohio/RR.pdf") #[c(2,8,15,21)]
spplot(Ohio.gen1, name2[c(2,8,15,21)], col.regions = brewer.pal(9, "Blues") , at=att,
       as.table=TRUE, main = "",colorkey = list(labels = list( labels = att,        
                                                               at = att)))
dev.off()

#at=c(0,0.2,0.4,0.6,0.8,0.9,1.0,2.6),

#P.RR
pdf("/home/edson/MCMC_SOU_3Ohio/p_RR.pdf") #[c(2,8,15,21)]
spplot(Ohio.gen1, name3[c(2,8,15,21)], col.regions = brewer.pal(7, "Greys") , at=c(0,0.2,0.4,0.6,0.8,1.0),
       as.table=TRUE, main = "",colorkey = list(labels = list( labels = c(0,0.2,0.4,0.6,0.8,1.0),        
                                                               at = c(0,0.2,0.4,0.6,0.8,1.0))))
dev.off()



pdf("/home/edson/MCMC_SOU_3Ohio/obsvpred.pdf")
plot(as.vector(x), as.vector(x.pred.mean), xlab="Observed", ylab="Predicted", 
     main="")
abline(lm(as.vector(x.pred.mean)~as.vector(x)))
#main=expression(paste(R^2, "= 0.99986"))
dev.off()

(cor(as.vector(x), as.vector(x.pred.mean)))^2


#Trend
dd1=read.table("/home/edson/MCMC_SOU_3Ohio/GDD_WOhioparamcb_beta.txt") #CHANGE ME
plot.ts(dd1[,175:182])
cc=dd1[2501:5000,181:268]
colnames(cc) = NN

aa = data$county.shp[1:n]
x.pred.mean1 = x.pred.mean2=x.pred.mean
cc.1 = cc.2 = cc
nname=numeric(n)
for (i in 1:n){
  x.pred.mean1[,aa[i]]=x.pred.mean2[,i]
  cc.1[,aa[i]]=cc.2[,i]
  nname[aa[i]]=colnames(x.pred.mean)[i]
}
colnames(cc.1)=nname

cc.2 = apply(cc.1,2,mean)
cc.2a = apply(cc.1, 2, function(x) quantile(x,c(0.025,0.975)))
cc.2l = cc.2a[1,]
cc.2u = cc.2a[2,]
sort(cc.2, decreasing=FALSE)
sort(cc.2l, decreasing=FALSE)

#Trend plots

pdf("/home/edson/MCMC_SOU_3Ohio/SMMRRRtrend.pdf")
plot(1968:1988,colMeans(RR.mean1.t), ylim=c(0.8,1.1), pch=1, xlab="", ylab="Risk estimate")
points(1968:1988,colMeans(SMR1.t),ylim=c(0.8,1), pch=2)
lines(1968:1988,1+(0.07*((1:21)/21))) #rescaled by a factor of 1
dev.off()


#Parameter extraction
#output[jj,]=c(2,0,taus2,tt[1],gamma,delta,beta2)
plot.ts(dd1[,1:10])
apply(dd1[2001:5000,c(3:5,50,92)],2,mean) #1,46,88 - Trumbull, Wayne and Butler
apply(dd1[2001:5000,c(3:5,50,92)],2,sd)
apply(dd1[2001:5000,c(3:5,50,92)],2,function(x) quantile(x,c(0.025,0.975)))

#b0, b1
apply(dd1[2501:5000,c(181,182)],2,mean)
apply(dd1[2501:5000,c(181,182)],2,sd)
apply(dd1[2501:5000,c(181,182)],2,function(x) quantile(x,c(0.025,0.975)))



#Other plots
RR.C1 = read.table("/home/edson/MCMC_SOU_3Ohio/CARBayesresults/OhioCARBayeRR_adaptive.txt", header=T)
RR.C1 = as.matrix(RR.C1)
fit.C1 = read.table("/home/edson/MCMC_SOU_3Ohio/CARBayesresults/OhioCARBAyefitted_adaptive.txt")
fit.C1 = data.matrix(fit.C1); fit.C1=as.vector(t(fit.C1))


RR.C2 = read.table("/home/edson/MCMC_SOU_3Ohio/CARBayesresults/OhioCARBayeRR_anova.txt", header=T)
RR.C2 = as.matrix(RR.C2)
fit.C2 = read.table("/home/edson/MCMC_SOU_3Ohio/CARBayesresults/OhioCARBAyefitted_anova.txt")
fit.C2 = data.matrix(fit.C2); fit.C2=as.vector(fit.C2)


RR.C3 = read.table("/home/edson/MCMC_SOU_3Ohio/CARBayesresults/OhioCARBayeRR_ar.txt", header=T)
RR.C3 = as.matrix(RR.C3)
fit.C3 = read.table("/home/edson/MCMC_SOU_3Ohio/CARBayesresults/OhioCARBAyefitted_ar.txt")
fit.C3 = data.matrix(fit.C3); fit.C3=as.vector(fit.C3)


RR.C4 = read.table("/home/edson/MCMC_SOU_3Ohio/CARBayesresults/OhioCARBayeRR_linear.txt", header=T)
RR.C4 = as.matrix(RR.C4)
fit.C4 = read.table("/home/edson/MCMC_SOU_3Ohio/CARBayesresults/OhioCARBAyefitted_linear.txt")
fit.C4 = data.matrix(fit.C4); fit.C4=as.vector(fit.C4)

RR.C5 = read.table("/home/edson/MCMC_SOU_3Ohio/CARBayesresults/OhioCARBayeRR_localised.txt", header=T)
RR.C5 = as.matrix(RR.C5)
fit.C5 = read.table("/home/edson/MCMC_SOU_3Ohio/CARBayesresults/OhioCARBAyefitted_localised.txt")
fit.C5 = data.matrix(fit.C5); fit.C5=as.vector(fit.C5)


#New vioplot
library(vioplot)
#pdf("Vioplots_Ohio.pdf")
pdf("/home/edson/MCMC_SOU_3Ohio/vio_obs_trend_all.pdf")
m <- rbind(c(1, 1, 1), c(2, 3, 4), c(5, 6, 7))
layout(m)
par(mar = c(3, 3, 1, 1), oma=c(1,1,0,0))
vioplot(as.vector(SMR1.t[,2:21]),as.vector(RR.mean1.t[,2:21]), as.vector(RR.C1[,2:21]), as.vector(RR.C2[,2:21]), as.vector(RR.C3[,2:21]), as.vector(RR.C4[,2:21]), as.vector(RR.C5[,2:21]),
        horizontal=FALSE, col="grey", ylim=c(0,3),names=c("SMR","Model 5", "CAR-1", "CAR-2", "CAR-3", "CAR-4", "CAR-5"))
mtext(text= "Risk estimate", side =2, line =2, cex=0.8)
#dev.off()  

#pdf("/home/edson/MCMC_SOU_3Ohio/obsvpred_all.pdf")
#par(mfrow=c(3,2))
plot(as.vector(x), as.vector(x.pred.mean), xlab="", ylab="", 
     main="Model 5")
mtext(text= "Predicted", side =2, line=2, cex=0.8)
abline(lm(as.vector(x.pred.mean)~as.vector(x)))
plot(as.vector(t(x)), fit.C1, xlab="", ylab="", 
     main="CAR-1")
abline(lm(fit.C1~as.vector(t(x))))
#main=expression(paste(R^2, "= 0.99986"))
plot(as.vector(t(x)), fit.C2, xlab="", ylab="", 
     main="CAR-2")
abline(lm(fit.C2~as.vector(t(x))))
#main=expression(paste(R^2, "= 0.99986"))
plot(as.vector(t(x)), fit.C3, xlab="", ylab="", 
     main="CAR-3")
mtext(text= "Predicted", side =2, line=2, cex=0.8)
mtext(text= "Observed", side =1, line=2, cex=0.8)
abline(lm(fit.C3~as.vector(t(x))))
#main=expression(paste(R^2, "= 0.99986"))
plot(as.vector(t(x)), fit.C4, xlab="", ylab="", 
     main="CAR-4")
mtext(text= "Observed", side =1, line=2, cex=0.8)
abline(lm(fit.C4~as.vector(t(x))))
#main=expression(paste(R^2, "= 0.99986"))
plot(as.vector(t(x)), fit.C5, xlab="", ylab="", 
     main="CAR-5")
mtext(text= "Observed", side =1, line=2, cex=0.8)
abline(lm(fit.C5~as.vector(t(x))))
#main=expression(paste(R^2, "= 0.99986"))
dev.off()


pdf("/home/edson/MCMC_SOU_3Ohio/Trend_obs_pred.pdf")
mat=matrix(c(1,1,2,3), 2, 2, byrow = TRUE)
#layout(mat, widths = rep.int(1, ncol(mat)), heights = rep.int(1, nrow(mat)))
layout(mat)
#plot1 - get data from Trend_plot.R
par(mar=c(4,3,1,2)) #bottom, left, top, and right
boxplot(cc[,1:88], ylab="", cex=0.8,names=1:88, cex.axis=1, col="grey")
mtext(text="Trend estimate", line=2, side=2, cex=0.8)
mtext(text="County", line=2, side=1, cex=0.8)
mtext(text="(a)", line=3, side=1, cex=0.8)

#plot2
par(mar=c(4,3,1,2)) #bottom, left, top, and right
plot(as.vector(x), as.vector(x.pred.mean), xlab="Observed \n (b)", ylab="")#main=expression(paste(R^2, "= 0.99986"))
mtext(text="Predicted", line=2, side=2, cex=0.8)
abline(lm(as.vector(x.pred.mean)~as.vector(x)))

#plot3
plot(as.vector(t(x)), fit.C, xlab="Observed \n (c)", ylab="")#main=expression(paste(R^2, "= 0.99768"))
mtext(text="Predicted", line=2, side=2, cex=0.8)
abline(lm(fit.C~as.vector(t(x))))
dev.off()



#Additional plots
par()
SMR.p = as.vector(SMR1.t)
RR.v = as.vector(RR.mean1.t)
plot(rep(0.2,1848),SMR.p, xlim=c(0,2), cex=1.5, xaxt="n", xlab="", ylab="", col="dodgerblue4")
mtext("Risk estimate", line=2,side=2,cex=1)
points(rep(0.8,1848),RR.v, cex=1.5, col="dodgerblue4")
for (i in 1:1848) segments(0.2, SMR.p[i], 0.8, RR.v[i], col= 'dodgerblue4')
axis(1,"SMR", at=0.2)
axis(1,"RR", at=0.8)

RR.v1 = as.vector(RR.C)
points(rep(1.2,1848),SMR.p, xlim=c(0,2), cex=1.5)
points(rep(1.8,1848),RR.v1, cex=1.5)
for (i in 1:1848) segments(1.2, SMR.p[i], 1.8, RR.v1[i], col= 'black')
axis(1,"SMR", at=1.2)
axis(1,"RR", at=1.8)


library(gridExtra)
pdf("/home/edson/MCMC_SOU_3Ohio/SMRRRPRplots.pdf") #name4[c(2,8,15,21)]
att=c(0,0.2,0.5,0.8,1.0,1.2,1.5,2.6)
grid.arrange(
  #SMR
spplot(Ohio.gen1, name4[c(2,8,15,21)], col.regions = brewer.pal(9, "Blues") , at=att,
            layout=c(4,1), as.table=TRUE, main = "",colorkey = list(labels = list( labels = att,        
                                                                 at = att)))
  ,
  
  #RR
spplot(Ohio.gen1, name2[c(2,8,15,21)], col.regions = brewer.pal(9, "Blues") , at=att,
            layout=c(4,1), as.table=TRUE, main = "",colorkey = list(labels = list( labels =att,        
                                                                 at = att)))
  
  ,
  #P.RR
spplot(Ohio.gen1, name3[c(2,8,15,21)], col.regions = brewer.pal(7, "Greys") , at=c(0,0.2,0.4,0.6,0.8,1.0),
            layout=c(4,1), as.table=TRUE, main = "", names.attr = c("1969","1975","1982", "1988"), colorkey = list(labels = list( labels = c(0,0.2,0.4,0.6,0.8,1.0),
                                                                 at = c(0,0.2,0.4,0.6,0.8,1.0))))        
)

dev.off()



plot(1968:1988,colMeans(RR.mean1.t), ylim=c(0.8,1.1), pch=1, xlab="", ylab="Risk estimate")
points(1968:1988,colMeans(SMR1.t),ylim=c(0.8,1), pch=2)
lines(1968:1988,1+(0.07*((1:21)/21))) #rescaled by a factor of 1

