Gtype_simulator<-function(MAF,n_Gtyp,n_iteration,plot){
  require(ggplot2)
  require(reshape2)
  require(gridExtra)
  itr<-0
  MAF_by_day<-cbind()
  genotype_pop2<-function(MAF,n_Gtyp){
    repeat{
      p<-1-MAF[1]
      q<-MAF[1]
      pop<-surv*200000 #this just simulates actual individuals. Current settings have genotype frequencies as output 
      #this number can be changed if necessary but very  low numbers will impair functionality
      ##############################################################################################
      #setting up genotype frequencies based on Hardy-Weinberg equillibrium and number of genotypes:
      ##############################################################################################
      if(n_Gtyp>2){  
        n_AA<-p^2*pop[1]
        n_AB<-2*p*q*pop[1]
        n_BB<-q^2*pop[1]
      }else{
        n_AA<-pop[1]/2
        n_AB<-pop[1]/2
        n_BB<-0
      }
      #creating output dataframe with headers:
      out<-data.frame("pop"=pop[1],"MAF"=q,
                      "Wbar"=NA,"Waa"=NA,"Wab"=NA,"Wbb"=NA,
                      "n_AA"=n_AA,"n_AB"=n_AB,"n_BB"=n_BB)
      #setting up starting (dummy) parameters:
      new_p<--1
      act_p<-1
      exp_AA<-exp_AB<-exp_BB<-pop[1]
      #calculating a random fitness value (W) for each genotype (aa,ab,bb)
      for(i in 1: (length(MAF)-1)){
        act_p<-1-MAF[i+1]
        p<-1-MAF[i]
        q<-MAF[i]
        ctr<-0
        Waa<-sample(0:100,1)/100
        Wab<-sample(0:100,1)/100
        Wbb<-sample(0:100,1)/100
        if(n_Gtyp<3){
          Wbb<-0 #because there are no bb genotypes if you only have 2 (aa/ab)
        }
        #calculating average fitness of the group:
        W_bar<-(p^2*Waa)+(2*p*q*Wab)+(q^2*Wbb)
        #based on the above (random) fitness, calculate expected change in p
        delta_p <-((p*q)/W_bar)*(p*(Waa-Wab)+q*(Wab-Wbb))
        new_p<- p + delta_p
        #calculate expected change in genotypes to next sampling point:
        p_AA<-(p^2*Waa)/W_bar
        p_AB<-(2*p*q*Wab)/W_bar
        p_BB<-(q^2*Wbb)/W_bar
        #calculate the number of individuals with each genotype:
        exp_AA<-p_AA * pop[i+1]
        exp_AB<-p_AB * pop[i+1]
        exp_BB<-p_BB * pop[i+1]
        #this step checks whether the *actual* number of individuals of each genotype fits with observed population size (survival)
        #**********************importantly, below is where you set the 'tolerance' of the program*****************************************
        #**********************i.e. how much error in counting of each genotype is allowable**********************************************
        #**********************default is set to 10% (i.e. 1.1 and 0.9) but can be modified if desired************************************
        #**********************some flexibility is necessary, but too much will be unhelpful**********************************************
        #**********************if the estimated change in number of each genotype is *possible* then it continues and W is re-estimated***
        #**********************if not, it tries again with new values of W to see if it can generate something that works*****************
        while((new_p > act_p*1.1 | new_p < act_p*0.9)|
              (exp_AA>n_AA*1.1|exp_AB>n_AB*1.1|exp_BB>n_BB*1.1)){
          Waa<-sample(1:100,1)/100
          Wab<-sample(1:100,1)/100
          Wbb<-sample(1:100,1)/100
          if(n_Gtyp<3){
            Wbb<-0
          }
          W_bar<-(p^2*Waa)+(2*p*q*Wab)+(q^2*Wbb)
          delta_p <-((p*q)/W_bar)*(p*(Waa-Wab)+q*(Wab-Wbb))
          new_p<- p + delta_p
          p_AA<-(p^2*Waa)/W_bar
          p_AB<-(2*p*q*Wab)/W_bar
          p_BB<-(q^2*Wbb)/W_bar
          exp_AA<-p_AA * pop[i+1]
          exp_AB<-p_AB * pop[i+1]
          exp_BB<-p_BB * pop[i+1]
          ctr<-ctr+1
          #if a feasible value for W cannot be found to satisfy the mortality and change in MAF within 10,000 tries, the program errors out,
          if(ctr>100){
            delta_p<-NA
            break
          }
        }
        if(!is.na(delta_p)){
          n_AA<-exp_AA
          n_AB<-exp_AB
          n_BB<-exp_BB
          
          new_pop<-sum(n_AA,n_AB,n_BB)
          out<-rbind(out,c(new_pop,round(1-new_p,3),
                           round(W_bar,3),round(Waa,3),round(Wab,3),round(Wbb,3),
                           round(n_AA),round(n_AB),round(n_BB)))
        }else{
          next
        }
      }
      if(nrow(out)==length(MAF)){
        return(out)
        break
      }
    }
  }
  for(a in 1:n_iteration){
    itr<-itr+1
    print(paste("running simulation number: ",itr))
    out<-genotype_pop2(MAF,n_Gtyp)
    if(itr==1){
      big_out<-out
    }else{
      big_out<-rbind(big_out,out)
    }
  }
  #assign days and percentage of each genotypeL
  big_out$Day<-s_days
  big_out$p_AA<-big_out$n_AA/(big_out$n_AA+big_out$n_AB+big_out$n_BB)
  big_out$p_AB<-big_out$n_AB/(big_out$n_AA+big_out$n_AB+big_out$n_BB)
  big_out$p_BB<-big_out$n_BB/(big_out$n_AA+big_out$n_AB+big_out$n_BB)
  if(plot=="Y"|plot=="YES"){
    #melt for graphing:
    delta_gtyp<-melt(big_out[,10:13],id.vars="Day")
    delta_fit<-melt(big_out[,c(10,4:6)],id.vars="Day")
    #adjust days for fitness (W) because each estimate is actually that of the sampling period prior:
    delta_fit<-na.omit(delta_fit)
    delta_fit$day<-as.numeric(lapply(delta_fit$Day,function(x){
      est_day<-grep(paste("^",x,"$",sep=""),s_days)
      act_day<-s_days[est_day-1]
      day<-as.numeric(act_day)
    }))
    #plot
    gtp<-ggplot(delta_gtyp,aes(Day,value*100))+
      geom_point(aes(color=variable),alpha=0.5,position = "jitter",na.rm = TRUE)+
      stat_smooth(aes(color=variable),method="loess",span=1)+
      scale_color_manual("Genotype",values=c("red2","green3","dodgerblue"))+
      scale_fill_manual("Genotype",values=c("red","green2","blue"))+
      geom_point(data=cbind.data.frame("af"=MAF,"day"=s_days),aes(x=day,y=af*100),
                 fill="purple",color="black",pch=23,size=5,alpha=0.5)+
      ylim(0,100)+
      xlim(-1,22)+
      labs(y="%",
           x="Day",
           title="Genotypes")+
      theme_minimal()
    fit<-ggplot(delta_fit,aes(day,value))+
      geom_point(aes(color=variable),alpha=0.5,position = "jitter",na.rm = TRUE)+
      stat_smooth(aes(color=variable),method="loess",span=1)+
      scale_color_manual("Fitness",values=c("red2","green3","dodgerblue"))+
      scale_fill_manual("FItness",values=c("red","green2","blue"))+
      ylim(0,1)+
      xlim(-1,22)+
      labs(y="Relative fitness (W)",
           x="Day",
           title="Fitness")+
      theme_minimal()
    require(gridExtra,)
    grid.arrange(gtp,fit)
    return(big_out)
  }
  else{
    return(big_out)
  }
}
