
## Abstract
Have you ever been on a cross country road trip, only to realise halfway through that you completely forgot to bring your Z table? Crash landed on an island, and needed normally distributed data to calibrate the raft you managed to build? Abducted by aliens that would only let you go on the condition that you give them data drawn from your species' favourite distribution? Studies have shown that 1 in 4 statisticians will experience a similar situation at some point in their life<sup>[1]</sup>. For this reason (and for fun), it's important to develop methods that can generate from something resembling a normal distribution in a way that's easy for a person to do without a calculator or computer. This article will explore different ways of generating normal data from simpler distributions that are easy to sample from, and transforming them using the Central Limit Theorem. Metrics will be needed to compare these different methods, so the Shapiro-Wilk normality test along with some basic sanity checks will be used.

## Introduction
Being able to sample from a normal distribution is a very common task in many statistical settings. The data can be used for simulations, as the proposal density for Markov Chain Monte Carlo, or to compare against other datasets. The word rnorm (R's random normal function) appears over 275,000 times in people's code on Github, meaning this kind of operation is quite common. Because of this, some very smart people have designed algorithms that can generate normally distributed samples incredibly quickly. For example, the ziggurat <sup>[2]</sup> and Box-Muller <sup>[5]</sup> methods can generate millions of samples per second. An issue with these types of algorithms is that they're designed for computers, so they make use of things that computers can do very well. They're not designed to be used by people. In terms of previous research, there doesn't seem to be much discussion of cost/efficiency tradeoffs, for cases where it's acceptable to lose some level of accuracy for a gain in speed. r-bloggers.com has a blog post that briefly explores generating a normal sample from 12 uniform samples<sup>[3]</sup>, which is related to the goal here, so it will be discussed and elaborated on.

## Standard Sampling Methods for Normal Distributions
The current state of the art appears to be the ziggurat method for sampling from a normal distribution. It is known to be faster than Box-Muller, although more complicated to implement. Computationally it's very fast, requiring only one random uniform number, a random index, a table lookup, a multiplication and a comparison. Although it's fast for a computer to do, it is difficult for a person to run. This is because the required lookup table is quite large, so it can't easily be memorised. The Julia programming language has an implementation of it <sup>[4]</sup> , and the constants used take up nearly 250 lines of code! Box-Muller is also frequently used (numpy uses it for their implementation). It doesn't require a lookup table and the equations are quite simple, but they involve logarithms, sines, cosines and square roots. If you have a calculator this is fine, but if you don't you would need to approximate each of these (maybe using a Taylor series), which can cause accuracy issues and skew the sampled distribution. Because people can't do math with 64 bits of precision or store gigabytes of data in working memory, it's important to base a sampling algorithm off of what we are good at: basic arithmetic and flipping coins.


## Testing Methodology
In order to distinguish between a good and a bad algorithm, we need some quantative way of telling how "close" a generated set of samples is to a normal distribution. One way of doing this is the Shapiro-Wilk test, which is a test statistic for this exact thing. More precisely, it tests against the null hypothesis that the data was generated from a normal distribution. Thus if the p value it produces is < 0.05, the null hypothesis can be rejected with 95% confidence, meaning that the data does not appear normal. So roughly, algorithms that give larger p values are more accurate. There are also a couple of sanity checks that can be performed, such as checking the the mean of the data is 0 and the variance is 1. Since what we're really trying to generate here is data that can pass as normal-ish, some less strict criteria will be used, such as checking if the data is symmetric, unimodal, and that the QQ plot in approximately linear. Doing calculations by hand can be slow, so the rough metric of "How long will this take" will also be taken into consideration.

The following methods rely on a few assumptions:
1. You can do basic arithmetic by hand
2. You can sample from a binomial distribution (you did bring your lucky quarter, right?)
3. You can sample from a uniform distribution, perhaps using a linear congruential generator (LCG)

## Binomial Distribution


First up is simply using coin flips. Assuming your coin isn't rigged, a Bernoulli trial using it will have mean 0.5 and variance 0.25. With this we can apply Central Limit Theorem to get that, asymptotically, (mean(n flips) - 0.5)/sqrt(0.25/n)) is normally distributed. The Central Limit Theorem is an incredibly beautiful result of statistics, but unfortunately it doesn't mean a whole lot at a face value since it doesn't say at what point the asymptotic behaviour kicks in. Using R, let's do some basic analysis on that. Starting with 5 flips, and sampling 1000 times:


```R
x = rep(0, 1000)
n = 5
for(i in 1:1000) {
  x[i] = (mean(rbinom(n,1,0.5)) - 0.5)/sqrt(0.25/n)
}
par(mfrow=c(1,2))
qqnorm(x, main=paste(n, "Flips"))
qqline(x)
hist(x, main=paste(n, "Flips"))
```


![png](Generating%20Approximately%20Normal%20Data%20by%20Hand_files/Generating%20Approximately%20Normal%20Data%20by%20Hand_5_0.png)


The QQ plot illustrates a very important shortcoming of this method: the binomial distribution is discrete, but a normal distribution is continuous! The mean can only have n+1 possible values (0-5 heads), which isn't an issue as n approaches infinity, but is very important if you don't want to tire out your thumb from flipping a coin constantly. Let's see how many flips it takes for the distribution to look acceptable, at least visually.


```R
x = rep(0, 1000)
n = c(25, 100, 1000)
par(mfrow=c(3,2))
# apply CLT to binomial samples 1000 times for 25, 100, and 1000 flips
for(i in 1:3) {
    for(j in 1:1000) {
        x[j] = (mean(rbinom(n[i],1,0.5)) - 0.5)/sqrt(0.25/n[i])
    }
    print(c(mean(x), var(x)))
    qqnorm(x, main=paste(n[i], "Flips"))
    qqline(x)
    hist(x, main=paste(n[i], "Flips"))
}
```

    [1] -0.0320000  0.9724284
    [1] 0.0474000 0.9885418
    [1] 0.01517893 0.97629790



![png](Generating%20Approximately%20Normal%20Data%20by%20Hand_files/Generating%20Approximately%20Normal%20Data%20by%20Hand_7_1.png)


All of the samples had approximately 0 mean and 1 variance, which is promising. Besides that, they don't look so great. 25 flips looks a bit choppy on the QQ plot, and the histogram seems to have much heavier tails, and overall not very normal. 100 flips looks better, but the ends of its QQ plot also shows that it has heavier tails than we'd like. 1000 flips looks quite good, but definitely fails from a practicality perspective. Based on this, flipping a coin to get a normal distribution seems to be a failed experiment. Failing generally isn't all that fun, but it's an important (and common) part of science. At least we figured this out now, as opposed to later down the line when you're desperately in need of random normal numbers and you try this first!

## Uniform Distribution
Next up is the sum of uniform variables. Sampling from a uniform distribution by hand isn't too difficult. You can use a linear congruential generator to get numbers that are pseudorandom enough for our use case. The majority (possibly all) LCGs in use have a modulus that's a power of 2 so they work better on computers, but I made one that's modulus 10000 for ease of computation (for the modulus all you need to do is grab the last 4 digits):


```R
lcg = function(seed) {
    return((61*seed + 7) %% 10000)
}
```

This satisfies the Hull-Dobell Theorem, so it will have a period of 10000 which should be sufficient for anything that can be done by hand. The next analysis will be very similar to what was done above, except using a uniform distribution. This will at the very least solve the issue where the binomial distribution was discrete, and will hopefully have quicker convergence! This will have some overlap in material with the blog post mentioned above, but it will be making use of this handcrafted LCG, and go deeper into the analysis. By the CLT, we have (sum(n samples from U(0,1)) - n*0.5)/sqrt(n/12) has an asymptotically normal distribution. When using n=12, this formula conveniently reduces to sum(n samples from U(0,1)) - 6. First off let's verify the result with n=12:


```R
x = rep(0, 1000)
seed = 42
par(mfrow=c(1,2))
# Grab 12 samples from the LCG and apply CLT 1000 times to get 1000 samples
for(i in 1:1000) {
  s = 0.0
  for(j in 1:12){
      seed = lcg(seed)
      s = s + seed/10000
  }
  x[i] = s - 6
}
qqnorm(x, main=paste(12, "Samples"))
qqline(x)
hist(x, main=paste(12, "Samples"))
shapiro.test(x)
library(moments)
print(paste("The mean of the data is ", mean(x)))
print(paste("The variance of the data is ", var(x)))
print(paste("The skew of the data is ", skewness(x)))
```


    
    	Shapiro-Wilk normality test
    
    data:  x
    W = 0.99942, p-value = 0.9933



    [1] "The mean of the data is  0.00460000000000002"
    [1] "The variance of the data is  1.02865793793794"
    [1] "The skew of the data is  -0.0246756170097312"



![png](Generating%20Approximately%20Normal%20Data%20by%20Hand_files/Generating%20Approximately%20Normal%20Data%20by%20Hand_12_2.png)


These results look good. The histogram looks normal, the data passes the Shapiro-Wilk normality test, and the mean, variance and skewness all seem to be quite close to what they should be. In terms of computation time, it actually isn't too bad. Drawing the 12 uniform samples only requires 12 multiplications and 12 additions, then add them together and subtract 6. You could even save some computation by saving 6 random samples from the previous run and drawing 6 new ones, which halves computation time at the risk of adding some correlation between the samples. While 12 seems like a good number to work with, it may be possible to use less samples at the loss of some accuracy.


```R
x = rep(0, 1000)
n = c(2,3,5,6)
par(mfrow=c(length(n),2))
seed = 42
data = as.table(matrix(rep(0, length(n)*5),ncol=5,byrow=TRUE))
colnames(data) <- c("# Samples", "Mean","Variance","Skewness", "SW P-value")
# Same as the above code, but for 2 3 5 and 6 samples
for(k in 1:length(n)){
    for(i in 1:1000) {
      s = 0.0
      for(j in 1:n[k]){
          seed = lcg(seed)
          s = s + seed/10000
      }
      x[i] = (s - n[k]*0.5)/sqrt(n[k]/12)
    }
    data[k, 1] = n[k]
    data[k, 2] = mean(x)
    data[k, 3] = var(x)
    data[k, 4] = skewness(x)
    data[k, 5] = shapiro.test(x)$p
    qqnorm(x, main=paste(n[k], "Samples"))
    qqline(x)
    hist(x, main=paste(n[k], "Samples"))
}
data
```


          # Samples          Mean      Variance      Skewness    SW P-value
    A  2.000000e+00  1.249240e-02  1.073739e+00 -5.931488e-02  1.144568e-05
    B  3.000000e+00 -1.670000e-02  1.057263e+00 -2.589752e-02  1.328073e-02
    C  5.000000e+00  4.260282e-03  1.008577e+00 -2.403329e-02  7.345421e-01
    D  6.000000e+00 -5.232590e-03  1.061236e+00 -1.026866e-01  1.826372e-01



![png](Generating%20Approximately%20Normal%20Data%20by%20Hand_files/Generating%20Approximately%20Normal%20Data%20by%20Hand_14_1.png)


All of these look approximately normal, in the sense that they're unimodal, symmetric about 0 with mean 0 and variance 1. Only 5 and 6 samples pass in this case, but if your current situation is okay with having slightly heavier tails, it's entirely possible to get away with using just 2 samples. One thing to be cautious about is that the random number generation can be sensitive to the LCG chosen. A key component is the length of the period. For numbers that are factors of 10000, this method will only generate 10000/n samples before repeating itself. If the number isn't a factor it can produce more before repeating, but there will be some overlap in the numbers used for new samples (for example one sample will use numbers 1-12 for a given seed, and another will use numbers 9996-9999 then 1-8). It's also possible that a particular dataset passing the Shapiro-Wilk test is highly sensitive to the initial seed. This is a good thing to check, or else all of the above results might not be representative. An easy way to verify this is the run the above code on a few hundred random seeds, and see what percent of the time each number of samples passes the test at a 0.05 threshold.


```R
x = rep(0, 1000)
n = c(2,3,4,5,6,12)
set.seed(42) # for reproducibility
# grab 500 random seeds, and them on each amount of samples
seeds = sample(1:10000, 500)
data = as.table(matrix(rep(0, length(n)*2),ncol=2,byrow=TRUE))
colnames(data) <- c("# Samples", "% of tests passed")

for(k in 1:length(n)){
    data[k, 1] = n[k]
    for(l in 1:length(seeds)){
        seed = seeds[l]
        for(i in 1:1000) {
          s = 0.0
          for(j in 1:n[k]){
              seed = lcg(seed)
              s = s + seed/10000
          }
          x[i] = (s - n[k]*0.5)/sqrt(n[k]/12)
        }
        if(shapiro.test(x)$p > 0.05) {
            data[k, 2] = data[k, 2] + 1
        }
    }
    data[k,2] = data[k,2] / length(seeds) * 100
}
data
```


      # Samples % of tests passed
    A       2.0               0.0
    B       3.0              29.6
    C       4.0              82.6
    D       5.0              97.0
    E       6.0              98.6
    F      12.0              99.8


These results are quite interesting. The fact that higher n is more likely to pass the Shapiro-Wilk test is in line with the Central Limit Theorem. What is more interesting is that this result doesn't seem to hold for all LCGs. For example, I accidentally discovered that using a = 81 instead and running the same analysis gives the following result:


```R
lcg2 = function(seed) {
    return((81*seed + 7) %% 10000)
}
x = rep(0, 1000)
n = c(2,3,4,5,6,12)
set.seed(42)
seeds = sample(1:10000, 500)
data = as.table(matrix(rep(0, length(n)*2),ncol=2,byrow=TRUE))
colnames(data) <- c("# Samples", "% of tests passed")

for(k in 1:length(n)){
    data[k, 1] = n[k]
    for(l in 1:length(seeds)){
        seed = seeds[l]
        for(i in 1:1000) {
          s = 0.0
          for(j in 1:n[k]){
              seed = lcg2(seed)
              s = s + seed/10000
          }
          x[i] = (s - n[k]*0.5)/sqrt(n[k]/12)
        }
        if(shapiro.test(x)$p > 0.05) {
            data[k, 2] = data[k, 2] + 1
        }
    }
    data[k,2] = data[k,2] / length(seeds) * 100
}
data
```


      # Samples % of tests passed
    A       2.0               0.0
    B       3.0              31.0
    C       4.0              92.2
    D       5.0              41.6
    E       6.0              99.4
    F      12.0              53.6


This is very much **not** in line with the Central Limit Theorem. The exact details explaining this result are most likely quite convoluted, but a possible theory is that since subsequent values from an LCG aren't actually independent, certain LCGs with certain sample sizes line up in such a way that averaging n samples results in some correlation that's negatively affecting the test. From a practical perspective, this isn't actually much of an issue. As long as your problem only requires samples that appear somewhat normal, you can memorise seeds that do seem to pass the Shapiro-Wilk test and use those.

Since 3 samples passes the test around 31% of the time, normal numbers can be generate quite cheaply. It takes a multiplication and addition to generate a number (the modulus and scaling division are both multiples of 10 so those barely count since it's just shifting numbers around), then 2 extra additions to add the 3 together. For scaling and shifting, 3 is a lucky number. Subraction by 1.5 is simple, and sqrt(3/12) = 1/2, so the scaling is simply multiplying by 2! So overall generating one random number takes 5 additions, 1 subtraction and 4 multiplications.

To perform this all you need to remember is the algorithm, the LCG, and a few seeds that can pass the Shapiro-Wilk test. Since the modulus is 10000, each seed can produce 3333 samples before the numbers used start overlapping. When picking seeds, you can really only pick 3 different ones, that aren't 3 apart from each other in the period or else it'll just repeat samples. Three such numbers are 6189, 7536 and 9703. 

This is most likely near the lower limit computation wise, at least for what a human can do without large precomputed tables or a calculator. It's possible that there's an LCG that can give passably normal results with 2 samples, which would reduce it to 3 additions, 1 subtraction and 3 multiplications, but that requires multiplying by sqrt(6), which will need an approximation so that can skew results.

With a calculator it's possible to do basic acceptance-rejection by sampling 2 uniform numbers, scaling the x value to perhaps the 99th quantile of the normal (around 2.3) and seeing if y < exp(-x^2/2). This is 3 multiplications, 2 additions and a few buttons on a calculator, but this only has an acceptance rate of around 43% so per valid sample this is actually more computation that 3 uniform variables, although more accurate. Without a calculator the exponential needs to be approximated with a Taylor series, and taking even 6 terms only approximates the normal well up to about 1.5 so being accurate would be too computationally intensive.

Overall we saw that flipping coins is an asymptotically correct way of sampling from a normal distribution, but impractical for a person to do. Sampling from a uniform distribution works quite well even when using less than 6 samples per normal number. For 3 or 12 samples the calculation doesn't involve inexact square roots, so they can be computed exactly. LCGs seem to be a relatively effective way to cheaply sample from a uniform distribution, although some of them, such as a = 81 and c = 7 from above have strange properties. Future research into this could try and find some precise way of describing this kind of behaviour. It would also be interesting to see if there's any LCG that can produce data that passes the Shapiro-Wilk test using only 2 uniform numbers, or if the pairwise correlation of LCGs prevents this from ever happening. Some further testing of the LCGs that do seem to work could also be useful, such as seeing if small subgroups appear normal, or if a large amount of samples is required for it to "average" out to a normal distribution.






<sup>[1]</sup> Citation needed

<sup>[2]</sup> http://www.jstatsoft.org/v05/i08/paper

<sup>[3]</sup> https://www.r-bloggers.com/clt-standard-normal-generator/

<sup>[4]</sup> https://github.com/JuliaLang/julia/blob/master/base/random.jl#L700

<sup>[5]</sup> http://projecteuclid.org/DPubS/Repository/1.0/Disseminate?view=body&id=pdf_1&handle=euclid.aoms/1177706645


```R

```
