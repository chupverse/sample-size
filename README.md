# RECQUIRED SAMPLE SIZES

<div style="text-align: justify">

This page is devoted to the calculation of the number of patients required for several studies such as an randomized clinical trial (RCT) for causal inference or a cohort for constructing or validating a predictive tool. R codes are proposed. We can click [here](https://poitiers-health-data.shinyapps.io/SampleSize/) to access to the related user-friendly calculators. Some reminders:

* **Superiority RCT:** Use to demonstrate that the experimental treatment is more effective than standard therapy.

* **Non-inferiority RCT:** Use to demonstrate that the experimental treatment is as effective as standard therapy.

* **Sequential RCT:** Intermediate analyses for early stopping the study.

* **Stepped wedge RCT:** Clusters are randomized sequentially to cross from control to experimental intervention.

## DESCRIBING A CHARACTERISTIC

<ul>
  <details>
  <summary>Mean</summary>
  
<br>
<em>
In order to describe an mean for an outcome with an expected standard deviation of 25 units with a total length of the 95% confidence interval equals to 10 units (5 units around the mean), the minimum sample size is 97 patients. (Chow, S.-C., Shao, J., Wang, H., & Lokhnygina, Y. (2017). Sample Size Calculations in Clinical Research (3rd ed.). Chapman and Hall/CRC.)
</em>

```r
sampleSize <- function(stDev, alpha, length)
{
Z <- qnorm(1-alpha/2)
return( (2 * Z * stDev / length)**2 )
}

sampleSize(stDev=25, alpha=0.05, length=10)

#> [1] 96.03647
```

**Input parameters:**
* stDev : expected standard deviation
* alpha : recquired type I error rate
* length : size of the (1-α)% confidence interval

  </details>	
</ul>

<ul>
  <details>
  <summary>Proportion</summary>
  
<br>
<em>
In order to describe an expected proportion of 35% with a total length of the 95% confidence interval equals to 10%, the minimum sample size is 350 patients. (Chow, S.-C., Shao, J., Wang, H., & Lokhnygina, Y. (2017). Sample Size Calculations in Clinical Research (3rd ed.). Chapman and Hall/CRC.)
</em>

```r
sampleSize <- function(p, alpha, length)
  {
  Z <- qnorm(1-alpha/2)
  return((((2*Z)**2)*(p*(1-p)))/(length**2))
  }

sampleSize(p=0.35, alpha=0.05, length=0.1)

#> [1] 349.5728
```

**Input parameters:**
* p : expected proportion of event
* alpha : recquired type I error rate
* length : total size of the (1-α)% confidence interval

  </details>
</ul>

## COMPARING TWO MEANS

### &nbsp;&nbsp;&nbsp;&nbsp;SUPERIORITY TRIALS

<ul>
  <details>
  <summary>No intermediate analysis : Individual randomization</summary>

<br>
<em>
Consider the following RCT with two parallel groups with a 1:1 randomization ratio. The expected mean is 66 units in patients in the experimental arm versus 72 units in the control arm. In order to demonstrate such a difference of 6 units, with a standard deviation of 23, a 5% two-sided type I error rate and a power of 80%, the minimum sample size per arm equals 231 (i.e., a total of 462 patients). (Chow, S.-C., Shao, J., Wang, H., & Lokhnygina, Y. (2017). Sample Size Calculations in Clinical Research (3rd ed.). Chapman and Hall/CRC.)
</em>

```r
library(epiR)

sample_mean <- function(mean1, mean0, sigma, power, r, sided.test, conf.level) {

  result <- epi.sscompc(treat = mean1, control = mean0,	sigma = sigma, n = NA, power = power, 
  		      	          r = r, sided.test = sided.test, conf.level = conf.level)
  		      	          
  if (result$n.treat<30 | result$n.control<30) {
    warning("--> At least one group size is < 30, normality assumption is questionnable and sample 
    size calculation may not be valid.")
  }

  return(result)
}

sample_mean(mean1 = 66, mean0 = 72, sigma = 23, power = 0.8, 
            r = 1, sided.test = 2, conf.level = 1-0.05)

#> $n.total
#> [1] 462

#> $n.treat
#> [1] 231

#> $n.control
#> [1] 231

#> $power
#> [1] 0.8

#> $delta
#> [1] 6
```

**Input parameters:**
* mean1: expected mean in the experimental arm
* mean0: expected mean in the control arm
* sigma: expected standard deviation in the two arms
* power: required power (1 minus type II error rate)
* r: randomization ratio (experimental:control)
* sided.test: one-sided test (1) or two-sided test (2) 
* conf.level: required confidence level (1 minus type I error rate)

  </details>
</ul>

<ul>
  <details>
  <summary>No intermediate analysis : Stepped wedge randomization</summary>
  
<br>
<em>
Consider the following stepped wedge RCT with 30 centers randomized in 30 sequences. The expected mean is 38 units in patients in the experimental arm versus 48 units in the control arm. In order to demonstrate such a difference of 10 units, with a standard deviation of 17 units, a 5% two-sided type I error rate and a power of 90%, the minimum sample size per arm equals 61 (i.e., a total of 122 patients) in case of individual randomization with a 1:1 ratio. According to our stepped wedge design and assuming an intraclass correlation coefficient of 0.05, we need to recruit 208 patients (104 in each arm). (Hemming K, Taljaard M. Sample size calculations for stepped wedge and cluster randomised trials: a unified approach. J Clin Epidemiol. 2016 Jan;69:137-46)
</em>

```r
library(epiR)

SampleSize_SW <- function(mean1, mean0, sigma, r, power, sided.test, conf.level,
center, sequence, icc) {
  
  SampSize_I <- epi.sscompc(treat = mean1, control = mean0, sigma = sigma, n = NA, 
                          r = r, power = power, sided.test = sided.test, conf.level = conf.level)
  
  if (SampSize_I$n.treat<30 | SampSize_I$n.control<30) {
    warning("--> At least one group size is < 30, in case of individual randomisation. Normality 
    assumption is questionnable and sample size calculation may not be valid.")
  }
  
  ni <- SampSize_I$n.total
  aa <- -2*center*(sequence - 1/sequence)*icc*(1+sequence/2)
  bb <- 3*ni*(1-icc)*icc*(1+sequence) - 2*center*(sequence -1/sequence)*(1-icc)
  cc <- 3*ni*(1-icc)*(1-icc) 
  m1 <- (-bb + sqrt(bb^2 - 4*aa*cc)) / (2*aa)
  m2 <- (-bb - sqrt(bb^2 - 4*aa*cc)) / (2*aa)
  m_sol <- max(m1,m2)
  Npat_center <- m_sol*(sequence+1) 
  N_tot_SW <- Npat_center*center 
  
  res <- list(SampSize_I$n.total,2*ceiling(N_tot_SW /2))
  names(res) <- c("n.indiv","n.SW")
  return(res)

}

SampleSize_SW(mean1 = 38, mean0 = 48, sigma = 17, r=1, power = 0.9, sided.test = 2,
              conf.level = 1-0.05, center = 30, sequence = 30, icc = 0.05)

# $n.indiv
# [1] 122
# 
# $n.SW
# [1] 208

```

**Input parameters:**
* mean1: expected mean in the experimental arm
* mean0: expected mean in the control arm
* sigma: expected standard deviation in the two arms
* r: individual randomization ratio (experimental:control)
* power: required power (1 minus type II error rate)
* sided.test: one-sided test (1) or two-sided test (2) 
* conf.level: required confidence level (1 minus type I error rate)
* center: number of centers
* sequence: number of sequences
* icc: expected intraclass correlation coefficient

	</details>
</ul>

<ul>
  <details>
	<summary>Sequential design</summary>
	
<br>
<em>
Consider the following RCT with two parallel groups with a 1:1 randomization ratio and 2 planned intermediate analyses for efficacy by using the O'Brien-Fleming method for considering the inflation of the type I error rate). The expected mean is 66 units in patients in the experimental arm versus 72 units in the control arm. In order to demonstrate such a difference of 6 units, with a standard deviation of 23, a 5% two-sided type I error rate and a power of 80%,  the final analysis should be carried out on 472 patients (236 patients per group). The first and second intermediate analyses would be performed on 158 and 316 patients respectively, i.e. 33% and 66% of the maximum number of included patients if their is no decision of stopping the study. (Demets DL, Lan KG. Interim analysis: The alpha spending function approach. Stat Med. 1994;13(13-14):1341–1352)
</em>

```r
library("rpact")
		
design <- getDesignGroupSequential(
               typeOfDesign = "OF", informationRates = c(1/3, 2/3, 1),
               alpha = 0.05, beta = 1-0.8, sided = 2)

designPlan <- getSampleSizeMeans(design, alternative = 6, stDev = 23,
                                 allocationRatioPlanned = 1)

summary(designPlan)

#> Stage                                          1       2       3 
#> Planned information rate                   33.3%   66.7%    100% 
#> Cumulative alpha spent                    0.0005  0.0143  0.0500 
#> Stage levels (two-sided)                  0.0005  0.0141  0.0451 
#> Efficacy boundary (z-value scale)          3.471   2.454   2.004 
#> Lower efficacy boundary (t)              -13.012  -6.405  -4.258 
#> Upper efficacy boundary (t)               13.012   6.405   4.258 
#> Cumulative power                          0.0329  0.4424  0.8000 
#> Number of subjects                         157.1   314.2   471.3 
#> Expected number of subjects under H1                       396.7 
#> Exit probability for efficacy (under H0)  0.0005  0.0138 
#> Exit probability for efficacy (under H1)  0.0329  0.4095 
```

**Input parameters:**
* typeOfDesign: type of design ("OF" for the O'Brien-Fleming method)
* informationRates: planned analyses defined as proportions of the maximum sample size
* alpha: recquired type I error rate
* beta: recquired type II error rate (1 minus power)
* sided: one-sided test (1), two-sided test (2)
* alternative: expected difference between the two arms
* stDev: expected standard deviation in the two arms
* allocationRatioPlanned: randomization ratio

</details>
</ul>

### &nbsp;&nbsp;&nbsp;&nbsp;NON-INFERIORITY TRIALS

<ul>
  <details>
  <summary>No intermediate analysis : Individual randomization</summary>

<br>
<em>
Consider the following RCT with two parallel groups with a 1:1 randomization ratio. The expected mean is 66 units in patients in the control arm and no difference compared to the experimental arm. Assuming an absolute non-inferiority margin of 7 points, a standard deviation of 23, the minimum sample size per arm equals 134 (i.e., a total of 268 patients) to achieve a 5% one-sided type I error rate and a power of 80%. (Chow, S.-C., Shao, J., Wang, H., & Lokhnygina, Y. (2017). Sample Size Calculations in Clinical Research (3rd ed.). Chapman and Hall/CRC.)
</em>

```r
library(epiR)

sample_mean <- function(mean0, sigma, delta, r, power, alpha) {

  result <- epi.ssninfc(treat = mean0, control = mean0,	sigma = sigma, delta = delta, n = NA,
                        r = r, power = power, alpha = alpha)
  		      	          
  if (result$n.treat<30 | result$n.control<30) {
  warning("--> At least one group size is < 30, normality assumption is questionnable and sample 
  size calculation may not be valid.")
  }

  return(result)
}

sample_mean(mean0 = 66, sigma = 23, delta = 7, r = 1, power = 0.8, alpha = 0.05)

#> $n.total
#> [1] 268

#> $n.treat
#> [1] 134

#> $n.control
#> [1] 134

#> $delta
#> [1] 7

#> $power
#> [1] 0.8
```
	
**Input parameters:**
* mean0: expected mean in both control and experimental arms
* sigma: expected standard deviation in the two arms
* delta: absolute non-inferiority margin
* r: randomization ratio (experimental:control)
* power: required power (1 minus type II error rate)
* alpha: required type I error rate


  </details>
</ul>

<ul>
  <details>
  <summary>No intermediate analysis : Stepped wedge randomization</summary>
  
<br>
<em>
Consider the following stepped wedge RCT with 30 centers randomized in 30 sequences. The expected mean is 48 units in patients in the control arm and no difference compared to the experimental arm. Assuming an absolute non-inferiority margin of 7 points, a standard deviation of 17, the minimum sample size per arm equals 102 (i.e., a total of 204 patients) to achieve a 5% one-sided type I error rate and a power of 90% in case of individual randomization with a 1:1 ratio. According to our stepped wedge design and assuming an intraclass correlation coefficient of 0.05, we need to recruit 372 patients (186 in each arm). (Hemming K, Taljaard M. Sample size calculations for stepped wedge and cluster randomised trials: a unified approach. J Clin Epidemiol. 2016 Jan;69:137-46)
</em>

```r
library(epiR)

SampleSize_SW <- function(mean0, sigma, delta, r, power, alpha, center, sequence, icc) {

  SampSize_I <- epi.ssninfc(treat = mean0, control = mean0, sigma = sigma, delta = delta, 
                            n = NA, r = r, power = power, alpha = alpha)
  
  if (SampSize_I$n.treat<30 | SampSize_I$n.control<30) {
    warning("--> At least one group size is < 30, in case of individual randomisation. Normality 
    assumption is questionnable and sample size calculation may not be valid.")
  }
  
  ni <- SampSize_I$n.total
  aa <- -2*center*(sequence - 1/sequence)*icc*(1+sequence/2)
  bb <- 3*ni*(1-icc)*icc*(1+sequence) - 2*center*(sequence -1/sequence)*(1-icc)
  cc <- 3*ni*(1-icc)*(1-icc) 
  m1 <- (-bb + sqrt(bb^2 - 4*aa*cc)) / (2*aa)
  m2 <- (-bb - sqrt(bb^2 - 4*aa*cc)) / (2*aa)
  m_sol <- max(m1,m2)
  Npat_center <- m_sol*(sequence+1) 
  N_tot_SW <- Npat_center*center 
  
  res <- list(SampSize_I$n.total,2*ceiling(N_tot_SW /2))
  names(res) <- c("n.indiv","n.SW")
  return(res)

}

SampleSize_SW(mean0 = 48, sigma = 17, delta = 7, r = 1, power = 0.9, alpha = 0.05,
              center = 30, sequence = 30, icc = 0.05)

# $n.indiv
# [1] 204
# 
# $n.SW
# [1] 372

```

**Input parameters:**
* mean0: expected mean in both control and experimental arms
* sigma: expected standard deviation in the two arms
* delta: absolute non-inferiority margin
* r: individual randomization ratio (experimental:control)
* power: required power (1 minus type II error rate)
* alpha: required confidence level (type I error rate)
* center: number of centers
* sequence: number of sequences
* icc: expected ntraclass correlation coefficient

  </details>
</ul>

<ul>
  <details>
  <summary>Sequential Design</summary>
  
<br>
<em>
This sample size is for a randomised controlled non-inferiority trial in two parallel groups experimental treatment versus control treatment with balanced randomisation (ratio 1 :1) for a continuous endpoint. Assuming an absolute non-inferiority margin of 7, with a standard deviation of 23, with a one-sided alpha risk of 5% and a power of 80%, the final analysis should be carried out on 276 patients(138 patients per group).Intermediate analyses would be performed on 92 and 184 patients respectively, i.e. 33%, 66% of the maximum number of included patients if their is no decision of stopping the study. (Demets DL, Lan KG. Interim analysis: The alpha spending function approach. Stat Med. 1994;13(13-14):1341–1352)
</em>

```r
library("rpact")
		
design <- getDesignGroupSequential(typeOfDesign = "OF", informationRates = c(1/3,2/3,1),
                                   alpha = 0.05, beta = 1-0.8, sided = 1)
                                   
designPlan <- getSampleSizeMeans(design, alternative = 0, stDev = 23,
                                 allocationRatioPlanned = 1, thetaH0 = -7)

summary(designPlan)

#> Stage                                          1       2       3 
#> Planned information rate                   33.3%   66.7%    100% 
#> Cumulative alpha spent                    0.0015  0.0187  0.0500 
#> Stage levels (one-sided)                  0.0015  0.0181  0.0437 
#> Efficacy boundary (z-value scale)          2.961   2.094   1.710 
#> Efficacy boundary (t)                      7.607   0.159  -2.246 
#> Cumulative power                          0.0660  0.4879  0.8000 
#> Number of subjects                          91.9   183.7   275.6 
#> Expected number of subjects under H1                       224.7 
#> Exit probability for efficacy (under H0)  0.0015  0.0172 
#> Exit probability for efficacy (under H1)  0.0660  0.4219 
```
**Input parameters:**
* typeOfDesign: type of design ("OF" for the O'Brien-Fleming method)
* informationRates: planned analyses defined as proportions of the maximum sample size
* alpha: recquired type I error rate
* beta: recquired type II error rate (1 minus power)
* sided: one-sided test (1)
* alternative: no difference between the two arms
* stDev: expected standard deviation in the two arms
* thetaH0 : equivalence limit
* allocationRatioPlanned: randomization ratio

  </details>
</ul>


## COMPARING TWO PROPORTIONS

### &nbsp;&nbsp;&nbsp;&nbsp;SUPERIORITY

<ul>
  <details>
    <summary>No intermediate analysis : Individual randomization</summary>
    
<br>
<em>
Consider the following RCT with two parallel groups with a 1:1 randomization ratio. The expected proportion of events is 35% in the experimental arm compared to 28% in the control arm. In order to demonstrate such a difference of 7%, with a two-sided type I error rate of 5% and a power of 80%, the minimum sample size per arm equals 691 (i.e., a total of 1,382 patients). (Chow, S.-C., Shao, J., Wang, H., & Lokhnygina, Y. (2017). Sample Size Calculations in Clinical Research (3rd ed.). Chapman and Hall/CRC.)
</em>

```r
library(epiR)

sample_proportion <- function(p1, p0, power, r, sided.test, conf.level) {

  result <- epi.sscohortc(irexp1 = p1, irexp0 = p0, n = NA, power = power,
                          r = r, sided.test = sided.test, conf.level = conf.level)
  
  if (result$n.exp1<30 | result$n.exp0<30) {
  warning("--> At least one group size is < 30, normality assumption is questionnable and sample 
  size calculation may not be valid.")
  }
  pmean <- (p1+p0)/2
  if (pmean*result$n.exp1 <5 | pmean*result$n.exp0 <5 
  | (1-pmean)*result$n.exp1 < 5 | (1-pmean)*result$n.exp0 <5) {
  warning("--> At least one theoretical effective is < 5, normality assumption is questionnable 
  and sample size calculation may not be valid.")
  }
  
  return(result)

}

sample_proportion(p1 = 0.35, p0 = 0.28, power = 0.80, 
                  r = 1, sided.test = 2, conf.level = 1-0.05)

#> $n.total
#> [1] 1382

#> $n.exp1
#> [1] 691

#> $n.exp0
#> [1] 691

#> $power
#> [1] 0.8

#> $irr
#> [1] 1.25

#> $or
#> [1] 1.384615
```
	
**Input parameters:**
*	p1: expected proportion in the experimental group
*	p0: expected proportion in the control group
*	power: required power (1 minus type II error rate)
* r: randomization ratio (experimental:control)
* sided.test: one-sided test (1), two-sided test (2)
* conf.level: recquired confidence level (1 minus type I error rate)

    </summary>
  </details>
</ul>


<ul>
  <details>
    <summary>No intermediate analysis : Stepped wedge randomization</summary>
    
<br>
<em>
Consider the following stepped wedge RCT with 15 centers randomized in 5 sequences. The expected proportion of events is 72% in the experimental arm compared to 62% in the control arm. In order to demonstrate such a difference of 10%, with a two-sided type I error rate of 5% and a power of 80%, the minimum sample size per arm equals 346 (i.e., a total of 692 patients) in case of individual randomization with a 1:1 ratio. According to our stepped wedge design and assuming an intraclass correlation coefficient of 0.01, we need to recruit 1,646 patients (823 in each arm). (Hemming K, Taljaard M. Sample size calculations for stepped wedge and cluster randomised trials: a unified approach. J Clin Epidemiol. 2016 Jan;69:137-46)
</em>

```r
library(epiR)

SampleSize_SW <- function(p1, p0, r, power, sided.test, conf.level, center, sequence, icc) {
  
  SampSize_I <- epi.sscohortc(irexp1 = p1, irexp0 = p0, n = NA, r = r,
                              power = power, sided.test = sided.test, conf.level = conf.level)
                              
  
  if (SampSize_I$n.exp1<30 | SampSize_I$n.exp0<30) {
  warning("--> At least one group size is < 30, normality assumption is questionnable and sample 
  size calculation may not be valid.")
  }
  pmean <- (p1+p0)/2
  if (pmean*SampSize_I$n.exp1 <5 | pmean*SampSize_I$n.exp0 <5 
  | (1-pmean)*SampSize_I$n.exp1 < 5 | (1-pmean)*SampSize_I$n.exp0 <5) {
  warning("--> At least one theoretical effective is < 5, normality assumption is questionnable 
  and sample size calculation may not be valid.")
  }

  ni <- SampSize_I$n.total
  aa <- -2*center*(sequence - 1/sequence)*icc*(1+sequence/2) 
  bb <- 3*ni*(1-icc)*icc*(1+sequence) - 2*center*(sequence -1/sequence)*(1-icc)
  cc <- 3*ni*(1-icc)*(1-icc)
  m1 <- (-bb + sqrt(bb^2 - 4*aa*cc)) / (2*aa)
  m2 <- (-bb - sqrt(bb^2 - 4*aa*cc)) / (2*aa)
  m_sol <- max(m1,m2) 
  Npat_center <- m_sol*(sequence+1) 
  N_tot_SW <- Npat_center*center 
  
  res <- list(SampSize_I$n.total,2*ceiling(N_tot_SW /2))
  names(res) <- c("n.indiv","n.SW")
  return(res)
}

SampleSize_SW(p1 = 0.72, p0 = 0.62, r = 1, power = 0.80, sided.test = 2, 
              conf.level = 1-0.05, center = 15, sequence = 5, icc = 0.01)

# $n.indiv
# [1] 692
# 
# $n.SW
# [1] 1646

```
	
**Input parameters:**
*	p1: expected proportion in the experimental group
*	p0: expected proportion in the control group
* r: individual randomization ratio (experimental:control)
*	power: required power (1 minus type II error rate)
* sided.test: one-sided test (1), two-sided test (2)
* conf.level: required confidence level (1 minus type I error rate)
* center: number of centers
* sequence: number of sequences
* icc: expected intraclass correlation coefficient

  </details>
</ul>

<ul>
  <details>
  <summary>Sequential design</summary>
  
<br>
<em>
Consider the following RCT with two parallel groups with a 1:1 randomization ratio and 2 planned intermediate analyses for efficacy by using the O'Brien-Fleming method for considering the inflation of the type I error rate. The expected proportion of event is 11% in patients in the experimental arm versus 15% units in the control arm. In order to demonstrate such a difference of 4%, with a 5% two-sided type I error rate and a power of 80%, the final analysis should be carried out on 2,256 patients (1,128 patients per group). The first and second intermediate analyses would be performed on 752 and 1,504 patients respectively, i.e. 33% and 66% of the maximum number of included patients if their is no decision of stopping the study. (Demets DL, Lan KG. Interim analysis: The alpha spending function approach. Stat Med. 1994;13(13-14):1341–1352)
</em>

```r
library("rpact")
		
design <- getDesignGroupSequential(typeOfDesign = "OF", 
                informationRates = c(1/3, 2/3, 1), alpha = 0.05,
                beta = 1-0.8, sided = 2)

designPlan <- getSampleSizeRates(design,  pi1 = 0.11, pi2 = 0.15,
                   allocationRatioPlanned = 1)

summary(designPlan)

#> Stage                                         1      2      3 
#> Planned information rate                  33.3%  66.7%   100% 
#> Cumulative alpha spent                   0.0005 0.0143 0.0500 
#> Stage levels (two-sided)                 0.0005 0.0141 0.0451 
#> Efficacy boundary (z-value scale)         3.471  2.454  2.004 
#> Lower efficacy boundary (t)              -0.079 -0.042 -0.029 
#> Upper efficacy boundary (t)               0.101  0.048  0.031 
#> Cumulative power                         0.0329 0.4424 0.8000 
#> Number of subjects                        751.8 1503.7 2255.5 
#> Expected number of subjects under H1                   1898.1 
#> Exit probability for efficacy (under H0) 0.0005 0.0138 
#> Exit probability for efficacy (under H1) 0.0329 0.4095 
```

**Input parameters:**
* typeOfDesign: type of design ("OF" for the O'Brien-Fleming method)
* informationRates: planned analyses defined as proportions of the maximum sample size
* alpha: required type I error rate
* beta: required type II error rate (1 minus power)
* sided: one-sided test (1), two-sided test (2)
* pi1: expected probability in the experimental group
* pi2: expected probability in the control group
* allocationRatioPlanned: randomization ratio (experimental/control)

  </details>
</ul>

### &nbsp;&nbsp;&nbsp;&nbsp;NON-INFERIORITY

<ul>
  <details>
    <summary>No intermediate analysis : Individual randomization</summary>
    
<br>
<em>
Consider the following RCT with two parallel groups with a 1:1 randomization ratio. The expected percentage of events is 35% in patients in the control arm and no difference compared to the experimental arm. Assuming an absolute non-inferiority margin of 5%,  the minimum sample size per arm equals 1,126 (i.e., a total of 2,252 patients) to achieve a 5% one-sided type I error rate and a power of 80%. (Chow, S.-C., Shao, J., Wang, H., & Lokhnygina, Y. (2017). Sample Size Calculations in Clinical Research (3rd ed.). Chapman and Hall/CRC.)
</em>

```r

sample_proportion <- function(p0, delta, r, power, alpha) {

  result <- epi.ssninfb(treat = p0, control = p0, delta = delta, n = NA, 
                          r = r, power = power, alpha = alpha)
  
  if (result$n.treat<30 | result$n.control<30) {
  warning("--> At least one group size is < 30, normality assumption is questionnable and sample size calculation may not be valid.")
  }
  pmean <- p0
  if (pmean*result$n.treat <5 | pmean*result$n.control <5 | (1-pmean)*result$n.treat < 5 | (1-pmean)*result$n.control <5) {
  warning("--> At least one theoretical effective is < 5, normality assumption is questionnable and sample 
  size calculation may not be valid.")
  }
  
  return(result)

}

sample_proportion(p0 = 0.35, delta = 0.05, r = 1, power = 0.8, alpha = 0.05)

#> $n.total
#> [1] 2252

#> $n.treat
#> [1] 1126

#> $n.control
#> [1] 1126

#> $delta
#> [1] 0.05

#> $power
#> [1] 0.8
```
	
**Parameters :**
* p0: expected proportion in both control and experimental arms
* delta: absolute non-inferiority margin
* r: randomization ratio (experimental:control)
* power: required power (1 minus type II error rate)
* alpha: required type I error rate


  </details>
</ul>

<ul>
  <details>
  <summary>No intermediate analysis : Stepped wedge randomization</summary>

<br>
<em>
Consider the following stepped wedge RCT with 15 centers randomized in 5 sequences. The expected proportion of events is 72% in patients in the control arm and no difference compared to the experimental arm. Assuming an absolute non-inferiority margin of 8%, the minimum sample size per arm equals 390 (i.e., a total of 780 patients) to achieve a one-sided type I error rate of 5% and a power of 80%, in case of individual randomization with a 1 :1 ratio. According to our stepped wedge design and assuming an intraclass correlation coefficient of 0.01, we need to recruit 1,890 patients (945 in each arm). (Hemming K, Taljaard M. Sample size calculations for stepped wedge and cluster randomised trials: a unified approach. J Clin Epidemiol. 2016 Jan;69:137-46)
</em>

```r
library(epiR)

SampleSize_SW <- function(p0, delta, r, power, alpha, center, sequence, icc) {
  
  SampSize_I <- epi.ssninfb(treat = p0, control = p0, delta = delta, n = NA, r = r,
                            power = power, alpha = alpha)
                              
  
  if (SampSize_I$n.treat<30 | SampSize_I$n.control<30) {
    warning("--> At least one group size is < 30, normality assumption is questionnable and sample 
    size calculation may not be valid.")
  }
  pmean <- p0
  if (pmean*SampSize_I$n.treat <5 | pmean*SampSize_I$n.control <5 
  | (1-pmean)*SampSize_I$n.treat < 5 | (1-pmean)*SampSize_I$n.control <5) {
    warning("--> At least one theoretical effective is < 5, normality assumption is questionnable 
    and sample size calculation may not be valid.")
  }

  ni <- SampSize_I$n.total
  aa <- -2*center*(sequence - 1/sequence)*icc*(1+sequence/2) 
  bb <- 3*ni*(1-icc)*icc*(1+sequence) - 2*center*(sequence -1/sequence)*(1-icc)
  cc <- 3*ni*(1-icc)*(1-icc)
  m1 <- (-bb + sqrt(bb^2 - 4*aa*cc)) / (2*aa)
  m2 <- (-bb - sqrt(bb^2 - 4*aa*cc)) / (2*aa)
  m_sol <- max(m1,m2) 
  Npat_center <- m_sol*(sequence+1) 
  N_tot_SW <- Npat_center*center

  res <- list(SampSize_I$n.total,2*ceiling(N_tot_SW /2))
  names(res) <- c("n.indiv","n.SW")
  return(res)
}

SampleSize_SW(p0=0.72, delta=0.08, r=1, power=0.8, alpha=0.05, center=15, sequence=5, icc=0.01)

# $n.indiv
# [1] 780
# 
# $n.SW
# [1] 1890
```

**Input parameters:**
* p0: expected proportion in both control and experimental arms
* delta: absolute non-inferiority margin
* r: individual randomization ratio (experimental:control)
* power: required power (1 minus type II error rate)
* alpha: required type I error rate
* center: number of centers
* sequence: number of sequences
* icc: expected intraclass correlation coefficient

  </details>
</ul>

<ul>
  <details>
  <summary>Sequential Design</summary>

<br>
<em>
This sample size is for a randomised controlled non-inferiority trial in two parallel groups experimental treatment versus control treatment with balanced randomisation (ratio 1 :1) for a binary endpoint. The expected percentage of events is 35% in patients in the control arm and no difference compared to the experimental arm. Assuming an absolute non-inferiority margin of 10%, with a one-sided alpha risk of 5% and a power of 80%, the final analysis should be carried out on 576 patients(288 patients per group).The two intermediate analyses would be performed on 192 and 384 patients respectively, i.e. 33%, 66% of the maximum number of included patients if their is no decision of stopping the study. (Demets DL, Lan KG. Interim analysis: The alpha spending function approach. Stat Med. 1994;13(13-14):1341–1352)
</em>

```r
library("rpact")
		
design <- getDesignGroupSequential(typeOfDesign = "OF", informationRates = c(1/3,2/3,1),
                                   alpha = 0.05, beta = 1-0.8, sided = 1)
                                   
designPlan <- getSampleSizeRates(design, pi1 = 0.35, pi2 = 0.35, thetaH0 = 0.10)

summary(designPlan)

#> Stage                                          1       2       3 
#> Planned information rate                   33.3%   66.7%    100% 
#> Cumulative alpha spent                    0.0015  0.0187  0.0500 
#> Stage levels (one-sided)                  0.0015  0.0181  0.0437 
#> Efficacy boundary (z-value scale)          2.961   2.094   1.710 
#> Efficacy boundary (t)                     -0.097  -0.002   0.032 
#> Cumulative power                          0.0660  0.4879  0.8000 
#> Number of subjects                         191.7   383.5   575.2 
#> Expected number of subjects under H1                       469.0 
#> Exit probability for efficacy (under H0)  0.0015  0.0172 
#> Exit probability for efficacy (under H1)  0.0660  0.4219 
```
**Input parameters:**
* typeOfDesign: type of design ("OF" for the O'Brien-Fleming method)
* informationRates: planned analyses defined as proportions of the maximum sample size
* alpha: recquired type I error rate
* beta: recquired type II error rate (1 minus power)
* sided: one-sided test (1)
* pi1 = pi2 : no difference between the two arms
* thetaH0 : equivalence limit
* allocationRatioPlanned: randomization ratio

  </details>
</ul>

## PREDICTING A PROPORTION

<ul>
  <details>
  <summary>Construction of predictive model</summary>
	
<br>
<em>
For developing a model/alghorithm based on 34 predictors as candidates with an expected R2 of at least 0.25 and an expected shrinkage of 0.9 (equation 11 in Riley et al. Statistics in Medicine. 2019;38:1276–1296), the minimal sample size is 1045.
</em>

```r
sampleSize <- function(predictors=34, R2=0.25, shrink=0.9)
 {  predictors/((shrink-1)*log(1-R2/shrink)) }

sampleSize()

#> [1] 1044.796
```

**Input parameters:**
* predictors : number of predictors as candidates
* R2 : expected R2
* shrink : expected shrinkage

  </details>
</ul>

<ul>
  <details>
  <summary>External validation</summary>
	
<br>
<em>
Consider O/E the ratio between the number of observed events versus expected ones. To achieve a precision defined as a length of the (1-α)% confidence interval of this ratio equals to 0.2, if the expected proportions is 50%, the required sample size is 386 (Riley et al. Minimum sample size for external validation of a clinical prediction model with a binary outcome. Statistics in Medicine. 2021;19:4230-4251).
</em>

```r
se <- function(width, alpha) # The standard error associated with the 1-alpha confidence interval
  {
  fun <- function(x) { exp( qnorm(1-alpha/2, mean=0, sd=1) * x ) - exp(-1* qnorm(1-alpha/2, mean=0, sd=1) * x ) - width } 
  return(uniroot(fun, lower = 0.001, upper = 100)$root)
  } 

size.calib <- function(p, width, alpha) # the minimum sample size to achieve this precision
  {   
  (1-p) / ((p * se(width=width, alpha=alpha)**2 ))
  }

size.calib(p=0.5, width=0.2, alpha=0.05)

#> [1] 385.4265
```
**Input parameters:**
* p: expected proportion of events
* width: size of the (1-α)% confidence interval
* alpha: type I error rate (α)

  </details>
</ul>
