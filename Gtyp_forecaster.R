gtp_4caster<-function(MAF,surv,n_Gtyp,pop){ 
  repeat{
    p<-1-MAF[1] 
    q<-MAF[1]
    pop<-pop * surv
    if(n_Gtyp>2){ #basic Hardy Weinberg assumptions
      n_AA<-p^2*pop[1]
      n_AB<-2*p*q*pop[1]
      n_BB<-q^2*pop[1]
    }else{     #if your known population is a cross between AA/AB parents:
      n_AA<-pop[1]/2
      n_AB<-pop[1]/2
      n_BB<-0
    }
    out<-data.frame("pop"=pop[1],"MAF"=q,
                    "Wbar"=NA,"Waa"=NA,"Wab"=NA,"Wbb"=NA,
                    "n_AA"=n_AA,"n_AB"=n_AB,"n_BB"=n_BB)
    new_q<--1 #set dummy variables to start
    act_q<-1
    exp_AA<-exp_AB<-exp_BB<-pop[1]
    
    for(i in 1: (length(MAF)-1)){
      act_q<-MAF[i+1] #this is the 'target' q at time t+1
      p<-1-MAF[i]     
      q<-MAF[i]
      ctr<-0
      while((new_q > act_q*1.1 | new_q < act_q*0.9)|
            (exp_AA>n_AA*1.1|exp_AB>n_AB*1.1|exp_BB>n_BB*1.1)){
        Waa<-sample(1:100,1)/100
        Wab<-sample(1:100,1)/100
        Wbb<-sample(1:100,1)/100
        if(n_Gtyp<3){
          Wbb<-0
        }
        W_bar<-(p^2*Waa)+(2*p*q*Wab)+(q^2*Wbb)             #calculate mean fitness
        delta_p <-((p*q)/W_bar)*(p*(Waa-Wab)+q*(Wab-Wbb))  #calculate change in p
        new_p<- p + delta_p                                #new simulated p at t+1
        new_q<- 1-new_p                                    #new simulated q at t+1  
        p_AA<-(p^2*Waa)/W_bar                              #new simulated genotype %
        p_AB<-(2*p*q*Wab)/W_bar
        p_BB<-(q^2*Wbb)/W_bar
        exp_AA<-p_AA * pop[i+1]                            #new simulated genotype #s
        exp_AB<-p_AB * pop[i+1]
        exp_BB<-p_BB * pop[i+1]
        ctr<-ctr+1
        if(ctr>10000){                                     #break after 10,000 tries
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
