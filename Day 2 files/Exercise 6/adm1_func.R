##################################################################
#####################Function to calculate adm2 level estimates
#####Save this as adm1_func.R
#ADM1 estimates and uncertainty (sd) 
#Note that pop is a population layer added to your pred_dat object
dd    <- 1:length(adm1_sp)
dd.un <- unique(spol1)
dmiss <- which(!dd%in%dd.un)

if (length(dmiss)>0) dd_num <- dd[-dmiss]
if (length(dmiss)==0) dd_num <- dd

dist_out <- matrix(0, length(dd_num), 5)
for (i in 1:length(dd_num)){
  if (length(which(spol1==dd_num[i]))==1){ 
    pop.ext <- pop[which(spol1==dd_num[i])] 
    ext <- as.vector(sapply(samp[which(spol1==dd_num[i]),], FUN=function(x) weighted.mean(x, w=pop.ext, na.rm=TRUE))) 
  }
  if (length(which(spol1==dd_num[i]))>1){  
    pop.ext <- pop[which(spol1==dd_num[i])]
    ext <- as.vector(apply(samp[which(spol1==dd_num[i]),], 2, FUN=function(x) weighted.mean(x, w=pop.ext, na.rm=TRUE)))
  }
  dist_out[i,] <- as.vector(c(mean(ext), sd(ext), quantile(ext, probs=c(0.025,0.5,0.975), na.rm=TRUE)))						
}

dist_out <- cbind(dd_num, dist_out)

