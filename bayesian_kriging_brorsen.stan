//Bayesian Kriging Brorsen
//Presented at StanCon 2020

data { 
  // Program uses a balanced panel - but can handle unbalanced
  int N; // number of obs 
  int M; // number of groups (years)
  int L; // number of locations
  
  vector[N] y; //the response variable
  matrix[L ,L] dist;
  int g[N];    // map obs to groups (years)
  int location[N]; // map obs to groups (locations)
  
  
}
parameters {
  
  
  real<lower=0,upper=100> a; //intercept
  real<lower=0,upper=15> rho; //diverges without the lower bound, but creates upward bias
  real<lower=1,upper=20> sill;  
  real<lower=0.1,upper=20> sigmae; 
  vector[L] b; //Spatial random effect plus mean
}
model {
  vector[N] local_mean; //intermediate result variable
  matrix[L,L] spatial_correlation;

  
  spatial_correlation=sill*exp(-rho*dist);

  b ~ multi_normal(rep_vector(a,L), spatial_correlation) ;
  
  for (j in 1:N) {
    local_mean[j]= b[location[j]];
  }
  
  
  //weakly informative priors, see section 6.9 in STAN user guide
  sill ~ gamma(2,.25); // expected value is the product of the first parameter times the inverse of the second
  sigmae~gamma(2,.2); // RSTAN uses inverse scale
  rho~normal(0,10);  
  a ~ normal(47,20);
  
  //model
  y ~ normal( local_mean,sigmae); //R uses standard deviation
  
}
