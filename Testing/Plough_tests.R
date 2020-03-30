####################################################################
### This script takes all biallelic markers from the dataset used by
### Plough and Hedgecock (2011)and compares simulated genotypes to #
### their estimated genotypic compositions #########################
####################################################################
rm(list=ls())
library(plyr)
library(reshape2)
library(ggplot2)
library(gridExtra)

#set up base survival and sampling days:
s_days<-c(0,2,6,10,16,22)
surv<-c(1,0.64,0.3,0.25,0.19,0.04)

#load all plough's 2 allele snps (with some larval sampling)
plough<-read.delim("/path/to/plough_mkrs.txt")
plough<-split(plough,plough$mkr)

#load the forecasting program
source("./path/to/Gtype_forecaster.R")

#check your inputs and set master lists of days and survival:
plough
surv
all_days<-s_days
all_surv<-surv
n<-3

#run the simulator on all markers in the list object:
sim_plough<-lapply(plough,function(x){
  n<<-3
  day_index<-sapply(x$Day,function(y){
    grep(paste("^",y,"$",sep=""),all_days)})
  s_days<<-all_days[day_index]
  surv<<-all_surv[day_index]
  if(max(x$p_BB)==0){
    n<<-2
    }
  Gtype_simulator(x$MAF,n,100,"N")
})

#check the outputs if desired
str(sim_plough,max.level = 1)
str(plough,max.level = 1)

#run another loop on the simulations to get plots:
ctr<-0
plots<-lapply(sim_plough,function(x){
  ctr<<-ctr+1
  m_mkr<-melt(x[,c(10:13)],id.vars="Day")
  m_ref<-melt(plough[ctr],id.vars=c("Day","MAF","mkr"))
  ggplot(m_mkr,aes(Day,value*100))+
    geom_point(aes(color=variable),alpha=0.3,position = "jitter")+
    stat_smooth(aes(color=variable),span=1)+
    scale_color_manual("Genotype",values=c("red2","green3","dodgerblue"),
                       labels=c("AA","AB","BB"))+
    geom_point(data=m_ref,aes(x=Day,y=value*100,fill=variable),color="black",pch=23,size=5)+
    scale_fill_manual("Genotype",values=c("red","green2","blue"),
                      labels=c("AA","AB","BB"))+
    ylim(0,100)+
    labs(y="% of pool",
         title = paste(as.character(unique(m_ref$mkr))),
         x=NULL)+
    guides(fill=FALSE,
           color=FALSE)+
    theme_minimal()
})

#plot the plots:
all_plot<-do.call("grid.arrange",c(plots,ncol=3))

#check the accuracy of the simulations:
acc_checker<-function(ac_freq,sim_freq,days){
  day_vec<-AA_vec<-AB_vec<-BB_vec<-c()
  for (a in days){
    pl_AA<-ac_freq$p_AA[ac_freq$Day==a]
    pl_AB<-ac_freq$p_AB[ac_freq$Day==a]
    pl_BB<-ac_freq$p_BB[ac_freq$Day==a]
    mkr<-as.character(unique(ac_freq$mkr))
    P.AA<-pnorm(pl_AA,mean=mean(sim_freq$p_AA[sim_freq$Day==a]),sd=sd(sim_freq$p_AA[sim_freq$Day==a]))
    P.AB<-pnorm(pl_AB,mean=mean(sim_freq$p_AB[sim_freq$Day==a]),sd=sd(sim_freq$p_AB[sim_freq$Day==a]))
    P.BB<-pnorm(pl_BB,mean=mean(sim_freq$p_BB[sim_freq$Day==a]),sd=sd(sim_freq$p_BB[sim_freq$Day==a]))
    day_vec<-append(day_vec,a)
    AA_vec<-append(AA_vec,round(P.AA,3))
    AB_vec<-append(AB_vec,round(P.AB,3))
    BB_vec<-append(BB_vec,round(P.BB,3))
  }
  return(cbind("day"=day_vec,
               "AA"=AA_vec,
               "AB"=AB_vec,
               "BB"=BB_vec,
               "mkr"=mkr))
}

acc_plough<-mapply(function(x,y){acc_checker(x,y,x$Day)},
  plough,sim_plough,SIMPLIFY = FALSE)

#bind all the data to one df:
all_acc<-as.data.frame(do.call("rbind",acc_plough))
all_acc$day<-as.numeric(as.character(all_acc$day))
all_acc$AA<-as.numeric(as.character(all_acc$AA))
all_acc$AB<-as.numeric(as.character(all_acc$AB))
all_acc$BB<-as.numeric(as.character(all_acc$BB))

#melt for plotting:
m_acc<-melt(all_acc,id.vars=c("day","mkr"),variable.name = "Genotype",
            value.name= "p.value")

#remove 'perfect' predictions (BB genotypes when none are present in 2-genotype cases)
m_acc<-m_acc[m_acc$p.value<1,]

#plot
ggplot(m_acc,aes(mkr,p.value))+
  geom_boxplot(aes(fill=Genotype))+
  scale_fill_manual(values=c("red","green3","dodgerblue"))+
  labs(x="Reference Marker",
       y="left-tailed p-value")+
  geom_hline(yintercept = c(0.95,0.5,0.05),
             color=c("red","black","red"),lty=c(1,2,1))+
  
  theme_minimal()

#how many fall outside a 95% CI?
nrow(m_acc[m_acc$p.value>0.95|m_acc$p.value<0.05,])/nrow(m_acc)


#Total magnitute of difference in mean simulations and empirical estimates:
av_sim<-mapply(function(x,y){
  av_sim<-ddply(x,"Day",summarize,
        meanAA=mean(p_AA),
        meanAB=mean(p_AB),
        meanBB=mean(p_BB))
  cbind(av_sim,y)},
  sim_plough,plough,SIMPLIFY = FALSE)

#bind the separate calculations:
av_sim<-do.call("rbind",av_sim)
av_sim$d_AA<-av_sim$meanAA-av_sim$p_AA
av_sim$d_AB<-av_sim$meanAB-av_sim$p_AB
av_sim$d_BB<-av_sim$meanBB-av_sim$p_BB

#melt for plotting:
m_av<-melt(av_sim[,c(8:14)],id.vars=c("MAF","Day","mkr"),
           variable.name = "Genotype",value.name="diff")

#plot:
ggplot(m_av,aes(diff*100))+
  geom_density(aes(fill=Genotype),alpha=0.5)+
  facet_wrap(~mkr,ncol=3)+
  scale_fill_manual(values=c("red","green3","dodgerblue"),
                    labels=c("AA","AB","BB"))+
  xlab("Difference in genotype frequency (%)")+
  theme_minimal()

#how many mean estimates fall outside of 10% of the 'true' value?
nrow(m_av[abs(m_av$diff)>0.1,])/nrow(m_av)
