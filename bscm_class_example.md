BSCM with Proposition 99 Data Set
================
Morgan Bale
2023-03-08

- <a href="#simulations" id="toc-simulations">Simulations</a>
- <a href="#smiulation-1-code-from-web-appendix"
  id="toc-smiulation-1-code-from-web-appendix">Smiulation 1: Code from Web
  Appendix</a>
- <a href="#simulation-2-code-from-paper"
  id="toc-simulation-2-code-from-paper">Simulation 2: Code from paper</a>
- <a href="#abadie-proposition-99-study"
  id="toc-abadie-proposition-99-study">Abadie Proposition 99 Study</a>

Load libraries and set your working directory to source file location
(we will work within the directory where you downloaded the files from
Github).

# Simulations

Kim et al (2020) provide two different version of their horseshoe BSCM
model. First, they provide equations in their paper that I turn into
stan code called, `bscm_paper_ex.stan`. Second, the authors provide a
web appendix which has stan code for their horseshoe model in section
B.1. I refer to this code from their web appendix as
`bscm_code_ex.stan`. The authors do not explain why they provide a
different parameterization of the model in the web appendix than in the
paper, so we test both here with simulated data. The paper
parameterization gives many more divergent transitions than the
parameterization from the web appendix.

Create simulated data to check recovery of the model code (you could
also do this in stan by sampling the parameters from their prior
distributions). The values are chosen to be similar to those provided in
the paper.

``` r
N_train=N_test=40
p = 5
beta_0=5
mu <- c(15, 35, 10, 20, 30)

X_train <- matrix(NA, nrow = N_train, ncol = p)
  
X_test <- matrix(NA, nrow=N_test, ncol=p)  #control unit matrix in post treatment 
  
for(j in 1:p) {
  X_train[,j] <- rnorm(N_train, mean=mu[j], sd=10) #control unit matrix in pre treatment
  X_test[, j] <- rnorm(N_test, mean=mu[j], sd=10)
}

beta <- c(-.5, 2, 0, 0, 0) 
  
#model 
epsilon <- rnorm(N_train, mean=0, sd=1)
y_train <- beta_0 + X_train%*%beta + epsilon

synthetic_data <- list(N_train=N_train, N_test=N_test, p=p, y_train=as.vector(y_train), X_train=X_train, X_test=X_test, beta=beta, beta_0=beta_0)
```

# Smiulation 1: Code from Web Appendix

Run model from web appendix provided by Kim et al (2020) with the
synthetic data created above. This model code is copied directly as is
from the web appendix.

Check results for simulation with web appendix model by looking at the
trace plots and values. Compare the printed values to the true values
chosen above.

``` r
print(bscm_code_model)
```

    ## S4 class stanmodel 'bscm_code_ex' coded as follows:
    ## // BSCM Code from Gupta Web Index
    ## // Morgan Bale 2023 
    ## 
    ## 
    ## data{
    ##     int N_train; //Number of observations in the pre-treatment periods
    ##     int N_test; //Number of observations in the post-treatment periods
    ##     int p; //Number of control units
    ##     real y_train[N_train]; //Treated unit in the pre-treatment periods
    ##     matrix[N_train, p] X_train; //Control unit matrix in the pre-treatment
    ##     matrix[N_test, p] X_test; //Control unit matrix in the post-treatment
    ## }
    ## 
    ## parameters{
    ##     real beta_0; //Intercept
    ##     real<lower=0> sigma2; //Error term variance
    ##     vector[p] beta_raw; //Control unit weights (will be transformed)
    ##     //Hyperparameters prior
    ##     vector<lower=0, upper=pi()/2>[p] lambda_unif;
    ##     real<lower=0> tau; //Global shrinkage
    ## }
    ## 
    ## transformed parameters{
    ##     vector[p] beta; //Control unit weights
    ##     real<lower=0> sigma; //Error term sd
    ##     vector<lower=0>[p] lambda; //Local shrinkage
    ##     vector[N_train] X_beta; //Synthetic control unit prediction in the pre-treatment period
    ##     lambda = tau * tan(lambda_unif); // => lambda ~ cauchy(0, tau)
    ##     for(j in 1:p){
    ##       beta[j] = lambda[j] * beta_raw[j];
    ##     }
    ##     sigma = sqrt(sigma2);
    ##     X_beta = beta_0 + X_train*beta;
    ## } 
    ## 
    ## model{
    ##     //Pre-treatment estimation
    ##     beta_raw ~ normal(0, 1); //=> beta ~ normal(0, lambda^2)
    ##     tau ~ cauchy(0, sigma);
    ##     sigma ~ cauchy(0,10);
    ##     beta_0 ~ cauchy(0,10);
    ##     y_train ~ normal(X_beta, sigma);
    ## }
    ## 
    ## generated quantities{
    ##   //Post-treatment prediction & Log-likelihood
    ##   vector[N_train] y_fit; //Fitted synthetic control unit in the pre-treatment
    ##   vector[N_test] y_test; //Predicted synthetic control unit in the post-treatment
    ##   vector[N_train] log_lik; //Log-likelihood
    ##   y_fit = beta_0 + X_train * beta;
    ##   for(i in 1:N_test){
    ##     y_test[i] = normal_rng(beta_0 + X_test[i,] * beta, sigma);
    ##   }
    ##   for (t in 1:N_train) {
    ##     log_lik[t] = normal_lpdf(y_train[t] | y_fit[t], sigma);
    ##   }
    ## }
    ## 
    ## 

``` r
traceplot(synth_draws, pars="beta")
```

![](bscm_class_example_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
traceplot(synth_draws, pars="beta_0")
```

![](bscm_class_example_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
traceplot(synth_draws, pars="sigma")
```

![](bscm_class_example_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r
traceplot(synth_draws, pars="lambda")
```

![](bscm_class_example_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->

``` r
#mcmc_recover_hist(As.mcmc.list(synth_draws, pars="beta"), true=as.vector(t(synthetic_data$beta)))
#mcmc_recover_hist(As.mcmc.list(synth_draws, pars="beta_0"), true=as.vector(t(synthetic_data$beta_0)))

print(synth_draws, pars="beta_0")
```

    ## Inference for Stan model: bscm_code_ex.
    ## 4 chains, each with iter=2000; warmup=1000; thin=1; 
    ## post-warmup draws per chain=1000, total post-warmup draws=4000.
    ## 
    ##        mean se_mean   sd 2.5%  25%  50% 75% 97.5% n_eff Rhat
    ## beta_0 6.04    0.05 0.83 4.46 5.47 6.04 6.6  7.65   248 1.03
    ## 
    ## Samples were drawn using NUTS(diag_e) at Thu Mar 23 14:21:14 2023.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

``` r
print(synth_draws, pars="beta")
```

    ## Inference for Stan model: bscm_code_ex.
    ## 4 chains, each with iter=2000; warmup=1000; thin=1; 
    ## post-warmup draws per chain=1000, total post-warmup draws=4000.
    ## 
    ##          mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
    ## beta[1] -0.50       0 0.02 -0.54 -0.52 -0.51 -0.49 -0.47  2496 1.00
    ## beta[2]  2.00       0 0.01  1.97  1.99  2.00  2.01  2.03   931 1.00
    ## beta[3]  0.00       0 0.01 -0.03 -0.01  0.00  0.00  0.02   591 1.01
    ## beta[4] -0.01       0 0.01 -0.04 -0.02 -0.01  0.00  0.02   457 1.02
    ## beta[5] -0.02       0 0.02 -0.06 -0.03 -0.02 -0.01  0.01   253 1.04
    ## 
    ## Samples were drawn using NUTS(diag_e) at Thu Mar 23 14:21:14 2023.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

Check the synthetic control plot as well.

``` r
#put synthetic control data together
y_fit <- summary(synth_draws, pars="y_fit")
sc_pre <- tibble(y_fit[[1]][,1])
lower <- y_fit[[1]][,4]
upper <- y_fit[[1]][,8]
sc_pre <- sc_pre %>% bind_cols(lower, upper)
```

    ## New names:
    ## • `` -> `...2`
    ## • `` -> `...3`

``` r
sc_pre <- sc_pre %>% mutate(week=rep(1:N_train))
names(sc_pre) <- c("synthetic_control","lower", "upper", "week")

#synthetic control in post period 
y_test <- summary(synth_draws, pars="y_test")
sc_post <- tibble(y_test[[1]][,1])
lower <- y_test[[1]][,4]
upper <- y_test[[1]][,8]
sc_post <- sc_post %>% bind_cols(lower, upper)
```

    ## New names:
    ## • `` -> `...2`
    ## • `` -> `...3`

``` r
sc_post <- sc_post %>% mutate(week=rep((N_train+1):(N_train+N_test)))
names(sc_post) <- c("synthetic_control", "lower", "upper", "week")

sc_data <- sc_pre %>% bind_rows(sc_post)

#now we need to add y 
y_test <- y_train + 10
y_data <- c(y_train, y_test)

#combine all data for plot
results_data <- sc_data %>% bind_cols(y_data)
```

    ## New names:
    ## • `` -> `...5`

``` r
names(results_data)[5] <- "generated_Y"

results_data %>% ggplot(aes(x=week)) + geom_ribbon(aes(ymin=lower, ymax=upper), fill="gray80") + geom_line(aes(y=generated_Y), color="darkred") + geom_line(aes(y=synthetic_control), color="steelblue") +
  labs(x="Week", y="Values") + ggtitle("Synthetic Control (blue) vs Generated Y (red)") + geom_vline(xintercept=N_train)
```

![](bscm_class_example_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# Simulation 2: Code from paper

Test with model code described in paper (different parameterization than
web appendix)

Check results for simulation with web appendix model (very similar
results, but the web appendix parameterization has less divergent
transitions)

``` r
print(bscm_paper_model)
```

    ## S4 class stanmodel 'bscm_paper_ex' coded as follows:
    ## // This model comes from Gupta's paper
    ## // Morgan Bale
    ## // March 2023
    ## 
    ## // Data
    ## data{
    ##   int N_train; //Number of observations in the pre-treatment periods
    ##   int N_test; //Number of observations in the post-treatment periods
    ##   int p; //Number of control units
    ##   real y_train[N_train]; //Treated unit in the pre-treatment periods
    ##   matrix[N_train, p] X_train; //Control unit matrix in the pre-treatment
    ##   matrix[N_test, p] X_test; //Control unit matrix in the post-treatment
    ## }
    ## 
    ## // The parameters accepted by the model. 
    ## parameters{
    ##   real<lower=0> sigma2; //Error term variance
    ##   vector[p] beta; 
    ##   //Hyperparameters prior
    ##   real<lower=0> tau; //Global shrinkage
    ##   real beta_0; //intercept 
    ##   vector<lower=0>[p] lambda; //Local shrinkage
    ## }
    ## 
    ## transformed parameters{
    ##   real<lower=0> sigma; //Error term sd
    ##   vector<lower=0>[p] lambda2; 
    ##   vector[N_train] X_beta; //Synthetic control unit prediction in the pre-treatment period
    ##   sigma = sqrt(sigma2);
    ##   X_beta = beta_0 + X_train*beta;
    ##   lambda2 = lambda .* lambda; 
    ## }
    ## 
    ## // The model to be estimated. 
    ## model{
    ##   //Pre-treatment estimation
    ##   beta ~ normal(0, lambda2);
    ##   lambda ~ cauchy(0, tau); 
    ##   tau ~ cauchy(0, sigma);
    ##   sigma ~ cauchy(0,10);
    ##   beta_0 ~ cauchy(0,10);
    ##   y_train ~ normal(X_beta, sigma);
    ## }
    ## 
    ## generated quantities{
    ##   //Post-treatment prediction & Log-likelihood
    ##   vector[N_train] y_fit; //Fitted synthetic control unit in the pre-treatment
    ##   vector[N_test] y_test; //Predicted synthetic control unit in the post-treatment
    ##   vector[N_train] log_lik; //Log-likelihood
    ##   y_fit = beta_0 + X_train * beta;
    ## 
    ##   for(i in 1:N_test){
    ##   y_test[i] = normal_rng(beta_0 + X_test[i,] * beta, sigma);
    ##     }
    ## 
    ##   for (t in 1:N_train) {
    ##   log_lik[t] = normal_lpdf(y_train[t] | y_fit[t], sigma);
    ##     }
    ## }

``` r
traceplot(synth_draws2, pars="beta")
```

![](bscm_class_example_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
traceplot(synth_draws2, pars="beta_0")
```

![](bscm_class_example_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
traceplot(synth_draws2, pars="sigma")
```

![](bscm_class_example_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->

``` r
traceplot(synth_draws2, pars="lambda")
```

![](bscm_class_example_files/figure-gfm/unnamed-chunk-7-4.png)<!-- -->

``` r
print(synth_draws2, pars="beta_0")
```

    ## Inference for Stan model: bscm_paper_ex.
    ## 4 chains, each with iter=2000; warmup=1000; thin=1; 
    ## post-warmup draws per chain=1000, total post-warmup draws=4000.
    ## 
    ##        mean se_mean   sd 2.5%  25% 50%  75% 97.5% n_eff Rhat
    ## beta_0 5.82    0.02 0.66  4.5 5.49 5.8 6.12  7.33   852    1
    ## 
    ## Samples were drawn using NUTS(diag_e) at Thu Mar 23 14:21:30 2023.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

``` r
print(synth_draws2, pars="beta")
```

    ## Inference for Stan model: bscm_paper_ex.
    ## 4 chains, each with iter=2000; warmup=1000; thin=1; 
    ## post-warmup draws per chain=1000, total post-warmup draws=4000.
    ## 
    ##          mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
    ## beta[1] -0.50    0.00 0.01 -0.53 -0.51 -0.50 -0.49 -0.48    15 1.10
    ## beta[2]  2.00    0.01 0.01  1.97  1.99  2.01  2.02  2.03     8 1.16
    ## beta[3]  0.00    0.00 0.01 -0.03 -0.01  0.00  0.00  0.01  1040 1.02
    ## beta[4] -0.01    0.01 0.01 -0.03 -0.02 -0.01  0.00  0.01     4 1.41
    ## beta[5] -0.02    0.00 0.02 -0.05 -0.03 -0.02  0.00  0.01    33 1.07
    ## 
    ## Samples were drawn using NUTS(diag_e) at Thu Mar 23 14:21:30 2023.
    ## For each parameter, n_eff is a crude measure of effective sample size,
    ## and Rhat is the potential scale reduction factor on split chains (at 
    ## convergence, Rhat=1).

Check synthetic control plot

``` r
#put synthetic control data together
y_fit <- summary(synth_draws2, pars="y_fit")
sc_pre <- tibble(y_fit[[1]][,1])
lower <- y_fit[[1]][,4]
upper <- y_fit[[1]][,8]
sc_pre <- sc_pre %>% bind_cols(lower, upper)
```

    ## New names:
    ## • `` -> `...2`
    ## • `` -> `...3`

``` r
sc_pre <- sc_pre %>% mutate(week=rep(1:N_train))
names(sc_pre) <- c("synthetic_control","lower", "upper", "week")

#synthetic control in post period 
y_test <- summary(synth_draws2, pars="y_test")
sc_post <- tibble(y_test[[1]][,1])
lower <- y_test[[1]][,4]
upper <- y_test[[1]][,8]
sc_post <- sc_post %>% bind_cols(lower, upper)
```

    ## New names:
    ## • `` -> `...2`
    ## • `` -> `...3`

``` r
sc_post <- sc_post %>% mutate(week=rep((N_train+1):(N_train+N_test)))
names(sc_post) <- c("synthetic_control", "lower", "upper", "week")

sc_data <- sc_pre %>% bind_rows(sc_post)

#now we need to add y 
y_test <- y_train + 10
y_data <- c(y_train, y_test)

#combine all data for plot
results_data <- sc_data %>% bind_cols(y_data)
```

    ## New names:
    ## • `` -> `...5`

``` r
names(results_data)[5] <- "generated_Y"

results_data %>% ggplot(aes(x=week)) + geom_ribbon(aes(ymin=lower, ymax=upper), fill="gray80") + geom_line(aes(y=generated_Y), color="darkred") + geom_line(aes(y=synthetic_control), color="steelblue") +
  labs(x="Week", y="Values") + ggtitle("Synthetic Control (blue) vs Generated Y (red)") + geom_vline(xintercept=N_train)
```

![](bscm_class_example_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

# Abadie Proposition 99 Study

Now that we have tested both code parameterizations from the BSCM model,
we will use the web appendix code version since there are less divergent
transitions. In this example, we will use the smoking data set from
Abadie et al (2010) with the BSCM provided by Kim et al (2020). This
allows us to compare the BSCM to the standard synthetic control method.
Abadie et al (2010) provides background on the data set and Proposition
99.

Data exploration: treatment time is 1988 making the pre period 1970 to
1988 and the post period 1989 to 2000

``` r
data(smoking)

length(unique(smoking$state)) #CA and 38 control states
```

    ## [1] 39

``` r
unique(smoking$year) #1970 to 2000
```

    ##  [1] 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984
    ## [16] 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999
    ## [31] 2000

``` r
#our y variable is cigsale (annual per capita cigarette consumption)
summary(smoking$cigsale)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    40.7   100.9   116.3   118.9   130.5   296.2

Prep data for stan

``` r
control_smoking <- smoking %>% filter(state!="California")
#make consistent order with states
control_smoking <-control_smoking %>% arrange(year,state)

treat_smoking <- smoking %>% filter(state=="California")
N_train = 19
N_test = 12
p = length(unique(control_smoking$state))

#(N_train, p)
X_train <- matrix(NA, nrow=N_train, ncol=p)

pre_period <- c(1970:1988)
i=1
for(n in pre_period) {
  X_train[i,] <- control_smoking %>% filter(year==n) %>% pull(cigsale)
  i=i+1
}

#(N_test, p)
X_test <- matrix(NA, nrow=N_test, ncol=p)

post_period <- c(1989:2000)
i=1
for(n in post_period) {
  X_test[i,] <- control_smoking %>% filter(year==n) %>% pull(cigsale)
  i=i+1
}

#length N_train
y_train <- treat_smoking %>% filter(year %in% c(1970:1988)) %>% pull(cigsale)

cig_data <- list(N_train=N_train, N_test=N_test, p=p, X_train=X_train, X_test=X_test, y_train=y_train)
```

Data exploration: compare to Table 1 and Figure 1 in Abadie et al
(2010). The recreation of Figure 1 shows why we don’t have an
appropriate control to use DID or another method.

``` r
rest_of_us <- control_smoking %>% group_by(year) %>% dplyr::summarise(avg_cigsale=mean(cigsale))

ca_cigsale <- treat_smoking %>% dplyr::select(year, cigsale) 
names(ca_cigsale) <- c("year", "CA")

rest_of_us %>% left_join(ca_cigsale, by="year") %>% ggplot(aes(x=year)) + geom_line(aes(y=avg_cigsale), color="steelblue") + geom_line(aes(y=CA), color="darkred") + labs(x="Year", y="Cigarette Sales") + ggtitle("Rest of US (blue) vs CA (red)") + geom_vline(xintercept=1988)
```

![](bscm_class_example_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
C1988 <- treat_smoking %>% filter(year==1988) %>% dplyr::summarise(mean(cigsale))
C1980 <- treat_smoking %>% filter(year==1980) %>% dplyr::summarise(mean(cigsale))
C1975 <- treat_smoking %>% filter(year==1975) %>% dplyr::summarise(mean(cigsale))

us1988 <- control_smoking %>% filter(year==1988) %>% dplyr::summarise(mean(cigsale))
us1980 <- control_smoking %>% filter(year==1980) %>% dplyr::summarise(mean(cigsale))
us1975 <- control_smoking %>% filter(year==1975) %>% dplyr::summarise(mean(cigsale))

table1 <- data.frame(vars=c("Sales per capita 1988", "Sales per capita 1980", "Sales per capita 1975"), California=c(as.character(C1988), as.character(C1980), as.character(C1975)), Avg_Controls=c(as.character(us1988), as.character(us1980), as.character(us1975)))

table1
```

    ##                    vars       California     Avg_Controls
    ## 1 Sales per capita 1988 90.0999984741211 113.823683688515
    ## 2 Sales per capita 1980 120.199996948242 138.089473724365
    ## 3 Sales per capita 1975 127.099998474121 136.931578987523

Run web appendix parameterization with cigarette data: since we did not
see a large difference between the paper and code parameterization, we
will use the web appendix code since it had less divergent transitions
(likely the reasons the authors used it).

Check traceplots: the traceplots suggest we might need more iterations
for better convergence (this can be adjusted in the sampling function).
We can also call individual beta and lambda traceplots to see them more
closely.

``` r
traceplot(cig_draws, pars="beta")
```

![](bscm_class_example_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
traceplot(cig_draws, pars="beta_0")
```

![](bscm_class_example_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

``` r
traceplot(cig_draws, pars="sigma")
```

![](bscm_class_example_files/figure-gfm/unnamed-chunk-13-3.png)<!-- -->

``` r
traceplot(cig_draws, pars="lambda")
```

![](bscm_class_example_files/figure-gfm/unnamed-chunk-13-4.png)<!-- -->

Check synthetic control plot (note the synthetic control is created in
generated quantities in the stan model). The Bayesian version allows us
to do the posterior shading to better understand our uncertainty (one of
the benefits of BSCM is that is provides a method for inference that
standard SC does not). This plot is recreating Figure 2 from Abadie et
al (2010).

``` r
#put synthetic control data together
y_fit <- summary(cig_draws, pars="y_fit")
sc_pre <- tibble(y_fit[[1]][,1])
lower <- y_fit[[1]][,4]
upper <- y_fit[[1]][,8]
sc_pre <- sc_pre %>% bind_cols(lower, upper)
```

    ## New names:
    ## • `` -> `...2`
    ## • `` -> `...3`

``` r
sc_pre <- sc_pre %>% mutate(year=rep(1970:1988))
names(sc_pre) <- c("synthetic_control","lower", "upper", "year")

#synthetic control in post period 
y_test <- summary(cig_draws, pars="y_test")
sc_post <- tibble(y_test[[1]][,1])
lower <- y_test[[1]][,4]
upper <- y_test[[1]][,8]
sc_post <- sc_post %>% bind_cols(lower, upper)
```

    ## New names:
    ## • `` -> `...2`
    ## • `` -> `...3`

``` r
sc_post <- sc_post %>% mutate(year=rep(1989:2000))
names(sc_post) <- c("synthetic_control", "lower", "upper", "year")

sc_data <- sc_pre %>% bind_rows(sc_post)

#now we need to add y 
y_data <- treat_smoking %>% dplyr::select(year, cigsale)
names(y_data) <- c("year", "CA")

#combine all data for plot
results_data <- sc_data %>% left_join(y_data, by="year")

results_data %>% ggplot(aes(x=year)) + geom_ribbon(aes(ymin=lower, ymax=upper), fill="gray80") + geom_line(aes(y=CA), color="darkred") + geom_line(aes(y=synthetic_control), color="steelblue") +
  labs(x="Year", y="Cigarette Sales") + ggtitle("Synthetic Control (blue) vs CA (red)") + geom_vline(xintercept=1988)
```

![](bscm_class_example_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

In Abadie et al (2010) they estimate that the mean decrease in packs per
capita was 20 and ours is about 15. We also recreate Figure 3 from their
paper.

``` r
results_data %>% mutate(diff=synthetic_control-CA) %>% filter(year %in% c(1989:2000)) %>% summarise(mean(diff))
```

    ## # A tibble: 1 × 1
    ##   `mean(diff)`
    ##          <dbl>
    ## 1         14.5

``` r
results_data %>% mutate(diff=CA-synthetic_control) %>% ggplot(aes(x=year, y=diff)) + geom_line() + labs(x="Year", y="Difference between SC and CA") + geom_hline(yintercept=0) + geom_vline(xintercept=1988)
```

![](bscm_class_example_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->
