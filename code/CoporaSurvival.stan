//Stan code based on ovarian corpora brms model
//adapted to include corpora to age data transformation
//using posterior distribution from demographic model 
//in Rand et al. (in review) Estimates of blue whale growth and reproduction parameters from historical data. Marine Mammal Science.
//Zoe Rand
//Last Updated: June 4, 2026
data {
  int<lower=1> N;  // total number of observations
  array[N] int Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> Kc;  // number of population-level effects after centering
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  int prior_only;  // should the likelihood be ignored?
  real omega; //for corpora rate prior
  real xi; //for corpora rate prior
  real alpha; //for corpora rate prior
}
transformed data {
 
}
parameters {
  vector[Kc] b;  // regression coefficients
  real Intercept;  // temporary intercept for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  array[M_1] vector[N_1] z_1;  // standardized group-level effects
  real Corpora_rate; //rate of corpora production
}
transformed parameters {
  vector[N_1] r_1_1;  // actual group-level effects
  // prior contributions to the log posterior
  real lprior = 0;
  r_1_1 = (sd_1[1] * (z_1[1]));
  lprior += student_t_lpdf(Intercept | 3, -2.3, 2.5);
  lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  //prior for corpora rate--mean 2.761, sd 0.424, alpha = 2.25 
  //but reparameterized in way that Stan wants
  lprior +=skew_normal_lpdf(Corpora_rate| xi, omega, alpha);
  
  //transforming data based on corpora rate
  matrix[N, K] NewAges; //new ages based on corpora rate
  for(i in 1:N){
    NewAges[i, 1]= X[i, 1];
    NewAges[i, 2] = X[i, 2]*Corpora_rate;
  }
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(NewAges[, i]);
    Xc[, i - 1] = NewAges[, i] - means_X[i - 1];
  }
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
    }
    target += poisson_log_glm_lpmf(Y | Xc, mu, b);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(z_1[1]);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  //posterior predictions
  // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    array[N] int post_pred; //save posterior predictions
    mu += Intercept;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
    }
    //target += poisson_log_glm_lpmf(Y | Xc, mu, b);
      post_pred = poisson_log_rng((mu + Xc*b));
  
}
