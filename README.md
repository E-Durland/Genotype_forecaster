# Genotype Forecaster

### Background

One of the major limitations to pooled DNA sequencing (pooseq) is that
genotype information is inherently sacrificed when the DNA of multiple
individuals is combined. The ‘raw material’ or output of poolseq is
allele frequency data, commonly reported as the relative abundance of
the less frequent or ‘minor’ allele (MAF).

In some scenarios, however, we can use additional information to infer
genotype composition. One such scenario is when a population is created
as a product of a controlled factorial cross (every male mated to every
female). In this case, genotype frequency in a pool of individuals
reflects the mendelian ratios resulting from the alleles of the gametes
that created it. In the associated script
(<a href="https://github.com/E-Durland/Genotype_simulator" class="uri">https://github.com/E-Durland/Genotype_simulator</a>)
we walked through how this understanding can be used to estimate the
genotypic composition of pools such as this.

In other settings, such as temporal re-sampling of a fixed population
size, the change in allele frequencies similarly observes certain
logical boundaries to how they will change. For example, if we have a
population of **n = 100** baby birds at **T = 0** with the genotypes:

``` r
birds
```

    ## AA AB BB 
    ## 75 15 10

the minor allele frequency is calculated as:

``` r
MAF <- (birds[[3]]*2 + birds[[2]]) / (sum(birds)*2)
MAF
```

    ## [1] 0.175

If we re-sample them at **T = 1**, and find **n = 30** with a **MAF =
0.30**, we lost 70% of our population but the ‘B’ allele rose in
frequency. In this case we know that we didn’t *‘gain’* any ‘BB’
genotypes to boost the MAF so we must have lost a lot of ‘AA’ birds
(since losing AB genotypes doesn’t change the relative frequency at
all). One possible combination is:

``` r
birds_t1
```

    ## AA AB BB 
    ## 20  2  8

``` r
MAF <- (birds_t1[[3]]*2 + birds_t1[[2]]) / (sum(birds_t1)*2)
MAF
```

    ## [1] 0.3

but there are several other combinations of genotypes that can sum to n
= 30 birds and MAF = 0.30. Instead of calculating each one by hand
(although it isn’t that hard), let’s just simulate all the combinations:

``` r
library(tidyr)
birds <-c ("AA"=75,"AB"=15,"BB"=10) #re-create starting population

AA <-AB <- BB <- c()  #make some empty objects
for(a in 1:50){   #create 50 simulations
  aa <- ab <- bb <- 0
    while(aa + ab + bb != 30 | #keep simulations to 30 birds 
          ((bb*2+ab)/((aa+ab+bb)*2)) != 0.3 | #keep sims for MAF = 0.30
          ( aa > birds[[1]] | ab > birds[[2]] | bb > birds[[3]] )){ #make sure each genotype doesn't exceed the max 
      aa<- sample(0:30,1) #randomly sample 1-30 birds for each genotype
      ab<- sample(0:30,1)
      bb<- sample(0:30,1)
    }
  AA <- append(AA, aa) #append 'good' simulations to objects
  AB <- append(AB, ab)
  BB <- append(BB, bb)
}
df<-cbind.data.frame(AA,AB,BB)  #make them a dataframe
df$MAF<- (df$BB*2+df$AB)/((df$AA+df$AB+df$BB)*2) #calculate MAF
df$ttl<-df$AA+df$AB+df$BB #calculate total # birds

birds_t1<-
  df[order(df$BB),]%>%
  unique.data.frame()
birds_t1
```

    ##    AA AB BB MAF ttl
    ## 7  14 14  2 0.3  30
    ## 8  15 12  3 0.3  30
    ## 4  16 10  4 0.3  30
    ## 1  17  8  5 0.3  30
    ## 2  18  6  6 0.3  30
    ## 3  19  4  7 0.3  30
    ## 10 20  2  8 0.3  30
    ## 15 21  0  9 0.3  30

With a single time step there are multiple ‘viable’ changes in genotype
combinations to satisfy empirical observations of mortality and changes
in MAF. With more time steps it becomes more difficult. If we re-sample
again at **T = 2** and find **n = 20** and **MAF = 0.20** what are the
possible scenarios? Taking just one:

``` r
ref_AA <- birds_t1$AA[3]
ref_AB <- birds_t1$AB[3]
ref_BB <- birds_t1$BB[3]
new_MAF <- 0.20 #this is the target MAF
surv <- 20 #make a max # of birds
AA <-AB <- BB <- c()
for(a in 1:10){   
  aa <- ab <- bb <- 0
  while(aa + ab + bb != surv | #keep simulations to **10** birds 
        ((bb*2+ab)/((aa+ab+bb)*2)) != new_MAF | #keep sims for MAF = target (0.20 in this case)
        ( aa > ref_AA | ab > ref_AB | bb > ref_BB )){ #make sure each genotype doesn't exceed the max 
    aa<- sample(0:surv,1) 
    ab<- sample(0:surv,1)
    bb<- sample(0:surv,1)
  }
  AA <- append(AA, aa) #append 'good' simulations to objects
  AB <- append(AB, ab)
  BB <- append(BB, bb)
}
df<-cbind.data.frame(AA,AB,BB)  #make them a dataframe
df$MAF<- (df$BB*2+df$AB)/((df$AA+df$AB+df$BB)*2) #calculate MAF
df$ttl<-df$AA+df$AB+df$BB #calculate total # birds
df[order(df$BB),]%>%
  unique.data.frame()
```

    ##    AA AB BB MAF ttl
    ## 10 12  8  0 0.2  20
    ## 2  13  6  1 0.2  20
    ## 9  14  4  2 0.2  20
    ## 1  15  2  3 0.2  20
    ## 4  16  0  4 0.2  20

We see that there are 5 combinations that satisfy these requirements. In
other scenarios (esp. with more individuals and less MAF change) it is
concievable that there are multiple possible genotype combinations that
can fit empirical observations. We can generalize these principals
through the Hardy Weinberg fitness equation:

> $\\overline{W}$ = *p*<sup>2</sup>*w*<sub>*A**A*</sub> +
> 2*p**q**w*<sub>*a**b*</sub> + *q*<sup>2</sup>*w*<sub>*B**B*</sub>

where p and q are the major and minor alleles, respectively,
*w*<sub>*x**x*</sub> = the fitness of each genotype and $\\overline{W}$
= mean fitness of the population. The change in the frequency of each
genotype is a function of it’s relative fitness (*w*<sub>*x**x*</sub>):

*F**r**e**q*<sub>*A**A*(*t*1)</sub> =
*F**r**e**q*<sub>*A**A*(*t*0)</sub> \* *w*<sub>*A**A*</sub>

and so on.

### The Function

In a system with sequential sampling of a fixed (upper-bounded)
population, we cannot know how genotypes change based on survival and
changes in MAF alone. With this data We can infer, however, what the
genotype changes *could not have been*. By randomly simulating fitness
for each genotype at each time point we can create a hypothesis of how
genotypes *might* be changing and reject them if they violate our
empirical observations (survival and MAF). Furthermore, we can simulate
these changes numerous times and create a range of hypothetical
trajectories from which to draw a mean and variance.

The core simulation functions have all been wrapped up in
[gtp\_4caster()](https://raw.githubusercontent.com/E-Durland/Genotype_forecaster/master/Gtyp_forecaster.R)

This script needs user inputs for:

-   MAF time series
-   Survival time series (cumulative)
-   Population starting size
-   \# genotypes at this locus

**gtp\_4caster** will run a single simulation through the given inputs.
For example:

``` r
source("/path/to/gtp_4caster.R")
MAF <- c(0.5,0.3,0.5)  #MAF changes from 50% to 30% to 50%
surv <- c(0.8,0.5,0.3) #Survival (cumulative) at 3 time points = 80%, 50%, 30%
n_Gtp <- 3             #This locus is expected to have 3 genotypes (AA/AB/BB)
pop <- 20000          #Starting population is 20k

#running the function:
gtp_4caster(MAF,surv,n_Gtp,pop)
```

    ##     pop   MAF  Wbar  Waa  Wab  Wbb n_AA n_AB n_BB
    ## 1 16000 0.500    NA   NA   NA   NA 4000 8000 4000
    ## 2 10000 0.320 0.618 0.93 0.75 0.04 3765 6073  162
    ## 3  6000 0.457 0.475 0.11 0.97 0.15  681 5148  171

*The NA’s in the first row of the output reflect that there is no
fitness value for the population at* **T = 0** *because mortality = 0*

This function can easily be put into a loop for multiple simulations:

``` r
# real data from an oyster larvae experiment:
surv<-c(1,0.64,0.30,0.25,0.19,0.04)  
pop<-200000
grad<-c(0.505102041,    0.387358401,    0.39731957, 0.423299638,    0.351726343,    0.440185334)
lin<-c(0.489795918, 0.229931981,    0.256375364,    0.280528109,    0.230838082,    0.285725147)
flip<-c(0.553030303,    0.407137554,    0.443218257,    0.42471454, 0.685093758,    0.375720974)

for(a in 1:100){
#  print(paste("Running simulation #",a))  #uncomment for verbose mode (to observe progress)
  sim<-gtp_4caster(flip,surv,3,pop)
  if(a == 1){
    out_df<-sim
  }else{
    out_df<-rbind(out_df,sim)
  }
}
out_df$Day<-rep(c(0,2,6,10,16,22),100) #append sampling times (days in this case)
out_df$p_AA<-out_df$n_AA/(out_df$n_AA+out_df$n_AB+out_df$n_BB) ##
out_df$p_AB<-out_df$n_AB/(out_df$n_AA+out_df$n_AB+out_df$n_BB) ##calculate percentages
out_df$p_BB<-out_df$n_BB/(out_df$n_AA+out_df$n_AB+out_df$n_BB) ##

m_pcent<-melt(out_df[,10:13],id.vars="Day")  #melt for plotting
ggplot(m_pcent,aes(Day,value*100))+
  geom_point(aes(color=variable),alpha=0.5,position = "jitter")+
  stat_smooth(aes(color=variable),span=1)+
  scale_color_manual("Genotype",values=c("red2","green3","dodgerblue"))+
  scale_y_continuous("% of population",limits=c(0,100))+
  theme_minimal()
```

![](Genotype_forecaster_files/figure-markdown_github/unnamed-chunk-9-1.png)
![](https://github.com/E-Durland/Genotype_forecaster/blob/master/Flip_ex.png)

For further examples and validation of this method, see the 'Testing' subdirectory
