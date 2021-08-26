#' #' A string that contains the R-Stan code for Bieri model 
#' with Gaussian likelihood
#' @param funcode: no parameter is required
#' @return The output is the corresponding R-Stan code
#' @examples
#' output <- fbieri_gauss();
#' @export
fbieri_gauss=function(){
content<-"functions { 
  vector fbieri(vector y, vector theta,
                              real[] x_r, int[] x_i){
 real x = y[1];
 real a = theta[1];
 real b = theta[2];
 real tmin = theta[3];
 real tmax = theta[4];
 real val;
 vector[1] y_x;
 val = a*(x)-a*tmin-pow(b,x-tmax);
 y_x[1] = val;
return y_x; }
  real my_normal_lpdf(real y, real mu, real sigma) {
  real seval;
    seval= -0.5*log(2*pi())-log(sigma)-(0.5*(y-mu)*(y-mu))/(sigma^2);
    return seval;
  }
}
data {
int<lower=0> N;
int<lower=0> NP; 
vector[NP] xpred;
vector[N] y;
vector[N] x;
real status[N];
}
transformed data{
vector[N] ey=exp(y);
real l1=min(x)-min(x)/20;
real l2=max(x)-max(x)/20;
vector[1] y_guess1;
vector[1] y_guess2;
real x_r[0];
int x_i[0];
//real st=1e-3;
y_guess1[1]= l1;
y_guess2[1]= l2;
}
parameters {
  real <lower=0, upper=1>  a;
real <lower=1>  b;
real <lower=0> tmin;
real <lower=tmin> tmax;
real <lower=0>  sigmasq;
}
transformed parameters {
real  sigma;
 vector[4] theta;
 theta[1]=a;
 theta[2]=b;
 theta[3]=tmin;
 theta[4]=tmax;
 sigma = sqrt(sigmasq);
}
model {
real ypred[N];
  real count=0;
for (i in 1:N)
{ypred[i] = a*(x[i])-a*tmin-pow(b,x[i]-tmax);
y[i] ~ normal(ypred[i], sigma);
count=x[i];
}
tmax ~ gamma(.1, .01); 
tmin ~ gamma(.1, .01); 
a ~ beta(1,1);
b ~ gamma(.2, .1);
sigmasq ~ inv_gamma(1e-3, 1e-3);
}
generated quantities { 
real log_lik[N];
real sscale2[NP];
real dev;
real tdmax;
real yeval[N];
real sh[N];
real ypred[NP];
real mu[NP];
vector[1] gmin;
vector[1] xmin;
vector[1] xmax;
vector[1] gmax;
gmin[1]=tmin;
gmax[1]=2*tmax;
for (i in 1:N) {
  sh[i] = a*(x[i])-a*tmin-pow(b,x[i]-tmax);  
log_lik[i] = normal_lpdf(y[i] |sh[i] , sigma);
yeval[i] = normal_rng(sh[i] , sigma);} 
for (j in 1:NP) {
sscale2[j] = a*(xpred[j])-a*tmin-pow(b,xpred[j]-tmax);  
if(is_nan(sscale2[j]))
{
  ypred[j]=0;
  sscale2[j]=0;
}else{
  ypred[j] = normal_rng(sscale2[j],sigma);
}
mu[j]=sscale2[j];}
dev = -2*sum(log_lik[]);
tdmax = (log(a)-log(log(b)))/log(b)+tmax;
xmin = algebra_solver(fbieri,gmin , theta, x_r, x_i, 1e-70, 1e+1, 1e+70); //for vbi
xmax = algebra_solver(fbieri, gmax , theta, x_r, x_i, 1e-70, 1e+1, 1e+70); //for vbi
} " 
return(content)
}
#' #' A string that contains the R-Stan code for Bieri model 
#' with Inverse Gamma likelihood. No zeros considered 
#' @param funcode: no parameter is required
#' @return The output is the corresponding R-Stan code
#' @examples
#' output <- fbieri_ig();
#' @export
fbieri_ig=function(){
  content<-"functions { 
  vector fbieri(vector y, vector theta,
                real[] x_r, int[] x_i){
    real x = y[1];
    real a = theta[1];
    real b = theta[2];
    real tmin = theta[3];
    real tmax = theta[4];
    real val;
    vector[1] y_x;
    val =a*(x-tmin)-pow(b,x)/pow(b,tmax);
    y_x[1] = val;
    return y_x;                              
  }
}
data {
  int<lower=0> N;
  int<lower=0> NP; 
  vector[NP] xpred;
  vector[N] y;
  vector[N] x;
}
transformed data{
  vector[N] ey=exp(y);
  real l1=min(x)-min(x)/2;
  real l2=max(x)+max(x)/2;
  vector[1] y_guess1;
  vector[1] y_guess2;
  real x_r[0];
  int x_i[0];
  y_guess1[1]= l1;
  y_guess2[1]= l2;
}
parameters {
  real <lower=0> tmax;
  real <lower=0, upper=1>  a;
  real <lower=1>  b;
  real <lower=0, upper=tmax> tmin;
  real <lower=2> shape;

}
transformed parameters{
  vector[4] theta;
  theta[1]=a;
  theta[2]=b;
  theta[3]=tmin;
  theta[4]=tmax;
}
model {
  vector[N]  mu;
  real count;
  real cntr=0; // center of the temperature
  real dst=0;
  for (i in 1:N)
  {
    mu[i] = a*(x[i])-a*tmin-pow(b,x[i]-tmax);
    if (!is_inf(mu[i]) && !is_nan(mu[i]) &&  mu[i]>0 )
      
          mu[i] = (shape-1)*mu[i];
          y[i] ~ inv_gamma(shape, mu[i]);
  }
  
  tmax ~ gamma(.1, .01); 
  tmin ~ gamma(.1, .01); 
  shape ~ gamma(.1, .01); 
  a ~ beta(1,1);
  b ~ gamma(.2, .1);
}
generated quantities { 
  real log_lik[N];
  real yeval[N];
  real sscale[N];
  real dev;
  real  <lower=tmin, upper=tmax> tdmax;
  vector[1] xmin;
  vector[1] xmax;
  real ypred[NP];
  real sscale2[NP];
  real mu[NP];
  
  real s2cale=0;
  xmin = algebra_solver(fbieri, y_guess1, theta, x_r, x_i, 1e-10, 1e+1, 1e+10);
  xmax = algebra_solver(fbieri, y_guess2, theta, x_r, x_i, 1e-10, 1e+1, 1e+10);
  for (i in 1:N) {
    sscale[i] = (a*(x[i])-a*tmin-pow(b,x[i]-tmax));
    if(!is_inf(sscale[i]) && !is_nan(sscale[i]) &&  sscale[i]>0 )
    {
      sscale[i] = (shape-1)*(sscale[i]);
      log_lik[i] = inv_gamma_lpdf(y[i] |shape , sscale[i]);
      yeval[i] = inv_gamma_rng(shape,sscale[i]);
    }else{
      log_lik[i] = 0;
      yeval[i] = 0;
    }
  }
  for (j in 1:NP) {
    sscale2[j] =(a*(xpred[j])-a*tmin-pow(b,xpred[j]-tmax)); 
    if(!is_inf(sscale2[j]) && !is_nan(sscale2[j]) &&  sscale2[j]>0 )
    {
      sscale2[j] = (shape-1)*(sscale2[j]);
      ypred[j] =log(inv_gamma_rng(shape,sscale2[j]));
      mu[j] = (a*(xpred[j])-a*tmin-pow(b,xpred[j]-tmax)); 
    }else{
      ypred[j]=0;
      sscale2[j] =0;
      mu[j] = 0;
    }}
  dev = -2*sum(log_lik[]);
  tdmax = (log(a)-log(log(b)))/log(b)+tmax;
} "
return(content)
}
#' #' A string that contains the R-Stan code for Bieri model 
#' with zero inflated Inverse Gamma likelihood. Excess zeros considered 
#' @param funcode: no parameter is required
#' @return The output is the corresponding R-Stan code
#' @examples
#' output <- fbieri_ziig();
#' @export
fbieri_ziig=function(){
  content<-"functions {
  vector fbieri(vector y, vector theta,
                              real[] x_r, int[] x_i){
 real x = y[1];
 real a = theta[1];
 real b = theta[2];
 real tmin = theta[3];
 real tmax = theta[4];
 real val;
 vector[1] y_x;
 val =a*(x-tmin)-pow(b,x)/pow(b,tmax);
 y_x[1] = val;
return y_x;
}
 real inv_logit_shift(real x, real r, real an){
 real val;
 real test=an*(x+r);
 vector[1] y_x;
  val=exp(- log_sum_exp(0,test));

return val;
}
  }
data {
  int<lower=0> N;
  int<lower=0> NP;
vector[NP] xpred;
vector[N] y;
vector[N] x;
vector[N] status;
}
transformed data{
real an=100;
real r=5e-3; // cause 0.05 is the median response and 0.1 is the max response.
real l1=min(x)-min(x)/2;
real l2=max(x)+max(x)/2;
vector[1] y_guess1;
vector[1] y_guess2;
real x_r[0];
int x_i[0];
y_guess1[1]= l1;
y_guess2[1]= l2;
}
parameters {
real <lower=0, upper=1>  a;
real <lower=1>  b;
real <lower=0> tmin; 
real <lower=tmin> tmax;
real <lower=2> shape;
}
transformed parameters{
  vector<lower=0, upper=1>[N] th;
 vector[4] theta;

 theta[1]=a;
 theta[2]=b;
 theta[3]=tmin;
 theta[4]=tmax;


for (i in 1:N) {
  th[i]=inv_logit_shift((a*(x[i])-a*tmin-pow(b,x[i]-tmax)),r,an);
}
}
model {
real log_other;
real log_zero;
vector[N]  ypred;

for (i in 1:N)
{
  if   (status[i]==0)   {
    log_zero=log(th[i]);
        ypred[i]=0;
       target += log_zero;
} else{
    ypred[i] = a*(x[i])-a*tmin-pow(b,x[i]-tmax);
  log_other=log(1-th[i]);
      target += log_other +
     inv_gamma_lpdf(y[i] | shape , (shape-1)* ypred[i]); //else ypred[i]=0;
}
}

tmin ~ gamma(.1, .01);
tmax ~ gamma(.1, .01);
shape ~ gamma(.1, .01);
a ~ beta(1,1);
b ~ gamma(.2, .1);
an ~ gamma(.1, .001);
}
generated quantities {
real log_lik[N];
real yeval[N];
real scale[N];
real dev;
real tdmax;
vector[1] xmin;
vector[1] xmax;
real ypred[NP];
real sscale2[NP];
real th2[NP];
real mu[NP];

real s2cale=0;
xmin = algebra_solver(fbieri, y_guess1, theta, x_r, x_i, 1e-10, 1e+1, 1e+10);
xmax = algebra_solver(fbieri, y_guess2, theta, x_r, x_i, 1e-10, 1e+1, 1e+10);


for (i in 1:N) {
scale[i] = (a*(x[i])-a*tmin-pow(b,x[i]-tmax));
 if (status[i]==0)
    {
    if (!is_nan(th[i])) log_lik[i] = bernoulli_lpmf(1 | th[i]); else log_lik[i] = 0;
    scale[i] = 0;
    }
    else if(!is_nan(scale[i]) && (!is_nan(th[i])) &&  scale[i]>0 && status[i]!=0)
    {
    log_lik[i] = bernoulli_lpmf(0 | th[i])+inv_gamma_lpdf(y[i] | shape,(shape-1)*scale[i]);
    } else
    {
    scale[i] = 0;
    log_lik[i] =0;
    }
    if(bernoulli_rng(th[i])|| scale[i]==0 || is_nan(scale[i])) yeval[i] =0; else yeval[i]=inv_gamma_rng(shape , scale[i]);
  }


for (j in 1:NP) {
th2[j]=inv_logit_shift((a*(xpred[j])-a*tmin-pow(b,xpred[j]-tmax)),r,an);
sscale2[j] =(a*(xpred[j])-a*tmin-pow(b,xpred[j]-tmax));

if(bernoulli_rng(th2[j]) || sscale2[j]<=0 || is_nan(sscale2[j]))
{
  ypred[j]=0;
  sscale2[j] =0;
    mu[j] =0;
}else{
   ypred[j] =(inv_gamma_rng(shape,(shape-1)*sscale2[j]));
  mu[j] = sscale2[j];

  }}
dev = -2*sum(log_lik[]);
tdmax = (log(a)-log(log(b)))/log(b)+tmax;
}"
return(content)
}
#' #' A string that contains the R-Stan code for Briere model 
#' with Gaussian likelihood
#' @param funcode: no parameter is required
#' @return The output is the corresponding R-Stan code
#' @examples
#' output <- fbriere_gauss();
#' @export
fbriere_gauss=function(){
  content<-"data {
int<lower=0> N;
real y[N];
real x[N];
real status[N];
int<lower=0> NP; 
vector[NP] xpred;
}
parameters {
real <lower=0> tmin; 
real <lower=tmin> tmax;
real <lower=0> a;
real  <lower=0> sigmasq;
}
transformed parameters {
real <lower=0>  sigma;
sigma = sqrt(sigmasq);
}
model {
real ypred[N];
 real ag;
  real bg;
  real cntr=0;
  real count=0;

for (i in 1:N)
{
  if (x[i] > tmin && x[i] < tmax ){ ypred[i] = exp(-a + log(x[i])+log(x[i]-tmin)+0.5*log(tmax-x[i]));
  y[i] ~ normal(ypred[i], sigma);
      count=x[i];
} else ypred[i]=0;
}
    tmax ~ gamma(.01, .001);
  tmin ~ gamma(.01,.01);
    a ~ gamma(0.1, 0.01);
sigmasq ~ inv_gamma(1e-3, 1e-3); //change from 3 to 1
}
generated quantities { 
real log_lik[N];
real yeval[N];
real sscale2[NP];
real ypred[NP];
real mu[NP];
real dev;
real alpha;
real  <lower=0> tdmax;
alpha = exp(-a);
for (i in 1:N) {
if ( x[i] > tmin && x[i] < tmax )
{
log_lik[i] = normal_lpdf(y[i] |alpha*(x[i])*(x[i]-tmin)*sqrt(tmax-x[i]), sigma);
yeval[i] = normal_rng(alpha*(x[i])*(x[i]-tmin)*sqrt(tmax-x[i]),sigma);
}
else
{
log_lik[i] = 0;
yeval[i] = 0;
}
}
for (j in 1:NP) {
sscale2[j] = alpha*(xpred[j])*(xpred[j]-tmin)*sqrt(tmax-xpred[j]);  
if(is_nan(sscale2[j]))
{
  ypred[j]=0;
  sscale2[j] =0;
}else{
  ypred[j] = normal_rng(sscale2[j],sigma);
}
mu[j]=sscale2[j];
}
dev = -2*sum(log_lik[]);
tdmax = ((4*tmax+3*tmin)+sqrt(pow((4*tmax+3*tmin),2)-40*tmin*tmax))/10;
} " 
return(content)
}
#' #' A string that contains the R-Stan code for Briere model 
#' with Inverse Gamma likelihood. No zeros considered 
#' @param funcode: no parameter is required
#' @return The output is the corresponding R-Stan code
#' @examples
#' output <- fbriere_ig();
#' @export
fbriere_ig=function(){
  content<-"data {
int<lower=0> N;
real y[N];
real x[N];
real status[N];
int<lower=0> NP; 
vector[NP] xpred;
}
parameters {
real <lower=0> tmin; 
real <lower=tmin> tmax;
real <lower=0> a;
  real <lower=2> shape;
}
model {
real mu[N];

for (i in 1:N)
{
    mu[i] = exp(-a+ log(x[i]) + log(x[i]-tmin) + (0.5)*log(tmax-x[i]));
  if(!is_nan(mu[i]) &&  mu[i]>0) target +=  inv_gamma_lpdf(y[i]|shape , (shape-1)*mu[i]);
}
  tmax ~ gamma(.01, .01);
  tmin ~ gamma(.01,.01);
  a ~ gamma(0.1, 0.01);
  shape ~ gamma(.01,.01);
}
 generated quantities { 
  real log_lik[N];
  real yeval[N];
  real dev;
  real tdmax;
  real scale[N];
  real alpha;
  int  count;
  real sscale2[NP];
  real ypred[NP];
  real mu[NP];
  count=0;
  alpha= exp(-a);
  for (i in 1:N) 
  {

              scale[i] = exp(-a)*x[i]*(x[i]-tmin)*sqrt(tmax-x[i]);
       if (!is_nan(scale[i]) &&  scale[i]>0 )//if (x[i] > tmin && x[i] < tmax )
    {
    scale[i]=(shape-1)*scale[i];
    log_lik[i] = inv_gamma_lpdf(y[i]|shape,scale[i]);
    yeval[i] = inv_gamma_rng(shape ,scale[i]);
    } else
    {
    scale[i] = 0;  
    log_lik[i] = 0;
    count +=1;
    yeval[i] = 0;
    }  
  }
  
  for (j in 1:NP) {
        sscale2[j] = exp(-a)*xpred[j]*(xpred[j]-tmin)*sqrt(tmax-xpred[j]);
  if ( !is_nan(sscale2[j]) &&  sscale2[j]>0 )  // if (xpred[j] > tmin && xpred[j] < tmax ) 
{
  
    ypred[j] = (inv_gamma_rng(shape,(shape-1)*sscale2[j]));
     mu[j] = sscale2[j];
}else{
  ypred[j]=0;
  sscale2[j] =0;
    mu[j] =0;
}}
  
  dev = -2*sum(log_lik[]);
tdmax = ((4*tmax+3*tmin)+sqrt(pow((4*tmax+3*tmin),2)-40*tmin*tmax))/10;
}  "
return(content)
}
#' #' A string that contains the R-Stan code for Briere model 
#' with zero inflated Inverse Gamma likelihood. Excess zeros considered 
#' @param funcode: no parameter is required
#' @return The output is the corresponding R-Stan code
#' @examples
#' output <- fbriere_ziig();
#' @export
fbriere_ziig=function(){
  content<-"functions{
  real[] fbriere(real tmin , real tmax , real a , real[] x ){
int N=num_elements(x);
real outl[N];
for(i in 1:N)
  {outl[i]= exp(-a + log(x[i]) + log(x[i] - tmin) + 0.5 * log(tmax - x[i]));
  }
  return outl;
  }
  real inv_logit_shift(real x, real r, real an){
 real val;
 real test=an*(x-r);
  val=exp(- log_sum_exp(0,test));

return val;
}
}
data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
  vector[N] status;
   int<lower=0> NP;
   vector[NP] xpred;
  }
  transformed data {
 real an=100;
real r=5e-3; //r=3.006212e-04; // cause 0.05 is the median response and 0.1 is the max response.
  }
  parameters {
  real <lower=0> tmin;
  real <lower=0> a;
  real <lower=tmin, upper=max(x)> tmax;
  real <lower=2> shape;
  }
  transformed parameters{
  vector<lower=0, upper=1>[N] th;
for (i in 1:N) {
  if ((x[i] > tmin) && (x[i] < tmax)) th[i]=inv_logit_shift(exp(-a+log(x[i])+log(x[i]-tmin)+(0.5)*log(tmax-x[i])),r,an); else th[i]=1;
}
}
  model {
  real ypred[N];
  real log_zero;
  real log_other;
  for (i in 1:N)
  {
      if (status[i]==0)
      {    log_zero=log(th[i]);
        target += log_zero; //only in case of zeros
            ypred[i]=0;
      }else
      {
    log_other=log(1-th[i]);
    ypred[i] = exp(-a+log(x[i])+log(x[i]-tmin)+(0.5)*log(tmax-x[i]));
    if ((!is_nan(ypred[i])) && (ypred[i]>0))  {
    target += log_other +
    inv_gamma_lpdf(y[i] | shape , (shape-1)* ypred[i]);
    } else  ypred[i]=0
    ;
    }
}
   tmin ~ gamma(.01,.01);
       a ~ normal(0,100);
        tmax ~ gamma(.01,.001);
    shape ~ gamma(.01, .001);
  }
  generated quantities {
  real log_lik[N];
  real yeval[N];
  real dev;
  real tdmax;
  real scale[N];
  real alpha;
  int  count;
   real mu[NP];
real ypred[NP];
real th2[NP];
  count=0;
  for (i in 1:N)
  {
        scale[i] = exp(-a+log(x[i])+log(x[i]-tmin)+(0.5)*log(tmax-x[i]));
if (status[i]==0)
    {
    if (!is_nan(th[i])) log_lik[i] = bernoulli_lpmf(1 | th[i]); else log_lik[i] = 0;
    scale[i] = 0;
    }
    else if(!is_nan(scale[i]) && (!is_nan(th[i]))&&  scale[i]>0 && status[i]!=0)
    {
    log_lik[i] = bernoulli_lpmf(0 | th[i])+inv_gamma_lpdf(y[i] | shape,(shape-1)*scale[i]);
    } else
    {
    scale[i] = 0;
    log_lik[i] =0;
    count +=1;
    }
    if(bernoulli_rng(th[i])|| scale[i]<=0 || is_nan(scale[i])) yeval[i] =0; else yeval[i]=inv_gamma_rng(shape , scale[i]);
  }

  for (j in 1:NP) {

    if ((xpred[j] > tmin) && (xpred[j] < tmax)) th2[j]=inv_logit_shift(exp(-a+ log(xpred[j])+log(xpred[j]-tmin)+(0.5)*log(tmax-xpred[j])),r,an); else th2[j]=1;
mu[j] = exp(-a+ log(xpred[j])+log(xpred[j]-tmin)+(0.5)*log(tmax-xpred[j]));
if(bernoulli_rng(th2[j]) || mu[j]<=0 || is_nan(mu[j]))
{
  ypred[j]=0;
    mu[j] =0;
}else{
   ypred[j] =(inv_gamma_rng(shape,(shape-1)*mu[j]));
  }}
  alpha= exp(-a);
  dev = -2*sum(log_lik[]);
tdmax = ((4*tmax+3*tmin)+sqrt(pow((4*tmax+3*tmin),2)-40*tmin*tmax))/10;
}"
  return(content)
}
#' #' A string that contains the R-Stan code for Analytis model 
#' with Gaussian likelihood
#' @param funcode: no parameter is required
#' @return The output is the corresponding R-Stan code
#' @examples
#' output <- fanalytis_gauss();
#' @export
fanalytis_gauss=function(){
  content<-"data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
  real status[N];
  int<lower=0> NP; 
vector[NP] xpred;
  }
  parameters {
      real  <lower=4> tmin;
 real  <lower=0> nn;
   real  <lower=0> mm;
  real <lower=35> tmax;
  real <lower=0> a;
  real <lower=0> sigmasq;
  }
   transformed parameters {
  real  sigma;
  sigma = sqrt(sigmasq);
  }
  model {
  real ypred[N];
  for (i in 1:N)
  {
 if (x[i] > tmin && x[i] < tmax)  
 {ypred[i] = exp(-a+ nn*log(x[i]-tmin)+mm*log(tmax-x[i])); 
    y[i] ~ normal(ypred[i], sigma);
}}
        tmin ~ gamma(.01,.01);
            nn ~ gamma(0.1, 0.1);
       mm ~ gamma(0.1, 0.1);
        tmax ~ gamma(.01,.001);
      a ~ gamma(0.1, 0.01);
sigmasq ~ inv_gamma(1e-3, 1e-3);}
  generated quantities { 
  real log_lik[N];
  real yeval[N];
  real dev;
  real tdmax;
  real scale[N];
  real alpha;
  int  count=0;
   real sscale2[NP];
   real mu[NP];
real ypred[NP];
  for (i in 1:N) 
  {
    if (x[i] > tmin && x[i] < tmax ) 
    {
    scale[i] = exp(-a+ nn*log(x[i]-tmin)+mm*log(tmax-x[i])); 
  log_lik[i] = normal_lpdf(y[i] | scale[i], sigma);
  yeval[i] = normal_rng(scale[i] , sigma);
    } else
    {
    scale[i] = 0;  
    log_lik[i] = 0;
    count +=1;
    yeval[i] = 0;
    }  
  }
    for (j in 1:NP) {
sscale2[j] = exp(-a+ nn*log(xpred[j]-tmin)+mm*log(tmax-xpred[j]));
if(is_nan(sscale2[j]))
{
  ypred[j]=0;
  sscale2[j]=0;
}else{
  ypred[j] = normal_rng(sscale2[j],sigma);
}
mu[j]=sscale2[j];}
  dev = -2*sum(log_lik[]);
  tdmax = (nn*tmax+mm*tmin)/(nn+mm);
    alpha= exp(-a);
} "
return(content)
}
#' #' A string that contains the R-Stan code for Analytis model 
#' with Inverse Gamma likelihood. No zeros considered 
#' @param funcode: no parameter is required
#' @return The output is the corresponding R-Stan code
#' @examples
#' output <- fanalytis_ig();
#' @export
fanalytis_ig=function(){
  content<-"data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
   int<lower=0> NP; 
vector[NP] xpred;
  }
  transformed data {
    vector[N] ey=exp(y);
  }
  parameters {
  real <lower=0> tmax;
  real  <lower=0> nn;
  real <lower=0> a;
  real  <lower=0> mm;
  real <lower=4, upper=tmax> tmin;
  real <lower=2> shape;
  

  
  }
  model {
  real ypred[N];
  real ag;
  real bg;
  real cntr=0;
  real count=0;
  for (i in 1:N)
  {
    ypred[i] = exp(-a+ nn*log(x[i]-tmin)+mm*log(tmax-x[i]));
    if(!is_nan(ypred[i]) &&  ypred[i]>0) target+= inv_gamma_lpdf(y[i]| shape , (shape-1)*(ypred[i]));
  }
        tmax ~ gamma(.01,.01);
    nn ~ gamma(0.1, 0.1);
       a ~ normal(0,100);
   tmin ~ gamma(.01,.01);
       mm ~ gamma(0.1, 0.1);
    shape ~ gamma(.01, .001); 
  }
  generated quantities { 
  real log_lik[N];
  real yeval[N];
  real dev;
  real tdmax;
  real scale[N];
  real alpha;
  int  count;
   real sscale2[NP];
   real mu[NP];
real ypred[NP];
  count=0;
  for (i in 1:N) 
  {
    scale[i] = exp(-a+nn*log(x[i]-tmin)+mm*log(tmax-x[i]));
    if (!is_nan(scale[i]) &&  scale[i]>0 )//if (x[i] > tmin && x[i] < tmax ) 
    {
    log_lik[i] = inv_gamma_lpdf(y[i] |shape , (shape-1)*scale[i]);
    yeval[i] = log(inv_gamma_rng(shape , (shape-1)*scale[i]));
    } else
    {
    scale[i] = 0;  
    log_lik[i] = 0;
    count +=1;
    yeval[i] = 0;
    }  
  }
   for (j in 1:NP) {
sscale2[j] = exp(-a+ nn*log(xpred[j]-tmin)+mm*log(tmax-xpred[j]));
  if ( !is_nan(sscale2[j]) &&  sscale2[j]>0 )  
{
    ypred[j] = (inv_gamma_rng(shape,(shape-1)*sscale2[j]));
     mu[j] = sscale2[j];
}else{
  ypred[j]=0;
    sscale2[j]=0;
    mu[j]=0;
  
}}
  
  
  alpha= exp(-a);
  dev = -2*sum(log_lik[]);
  tdmax = (nn*tmax+mm*tmin)/(nn+mm);
}  "
return(content)
}
  #' #' A string that contains the R-Stan code for Analytis model 
  #' with zero inflated Inverse Gamma likelihood. Excess zeros considered 
  #' @param funcode: no parameter is required
  #' @return The output is the corresponding R-Stan code
  #' @examples
  #' output <- fanalytis_ziig();
  #' @export
  fanalytis_ziig=function(){
    content<-"functions{
  real inv_logit_shift(real x, real r, real an){
 real val;
 real test=an*(x-r);
  val=exp(- log_sum_exp(0,test));

return val;
}
}
data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
  vector[N] status;
   int<lower=0> NP;
   vector[NP] xpred;
  }
  transformed data {
 real an=100;
real r=5e-3; //r=3.006212e-04; // cause 0.05 is the median response and 0.1 is the max response.
}
  parameters {
  real <lower=4> tmin;
  real  <lower=0> nn;
  real <lower=0> a;
  real  <lower=0> mm;
  real <lower=tmin> tmax;
  real <lower=2> shape;
  }
  transformed parameters{
  vector<lower=0, upper=1>[N] th;


for (i in 1:N) {

  if ((x[i] > tmin) && (x[i] < tmax)) th[i]=inv_logit_shift(exp(-a+nn*log(x[i]-tmin)+mm*log(tmax-x[i])),r,an); else th[i]=1;
}
}
  model {
  real ypred[N];
  real log_zero;
  real log_other;
  for (i in 1:N)
  {

      if (status[i]==0)
      {    log_zero=log(th[i]);
        target += log_zero; //only in case of zeros
            ypred[i]=0;
      }else
      {
    log_other=log(1-th[i]);
    ypred[i] = exp(-a+nn*log(x[i]-tmin)+mm*log(tmax-x[i]));
    if ((!is_nan(ypred[i])) && (ypred[i]>0))  {
    target += log_other +
    inv_gamma_lpdf(y[i] | shape , (shape-1)* ypred[i]);
    } else  ypred[i]=0;
    }
}
   tmin ~ gamma(.01,.01);
    nn ~ gamma(0.1, 0.1);
       a ~ normal(0,100);
        tmax ~ gamma(.01,.001);
       mm ~ gamma(0.1, 0.1);
    shape ~ gamma(.01, .001);
  }
  generated quantities {
  real log_lik[N];
  real yeval[N];
  real dev;
  real tdmax;
  real scale[N];
  real alpha;
  int  count;
   real mu[NP];
real ypred[NP];
real th2[NP];

  count=0;
  for (i in 1:N)
  {

        scale[i] = exp(-a+nn*log(x[i]-tmin)+mm*log(tmax-x[i]));
if (status[i]==0)
    {
    if (!is_nan(th[i])) log_lik[i] = bernoulli_lpmf(1 | th[i]); else log_lik[i] = 0;
    scale[i] = 0;
    }
    else if(!is_nan(scale[i]) && (!is_nan(th[i]))&&  scale[i]>0 && status[i]!=0)
    {
    log_lik[i] = bernoulli_lpmf(0 | th[i])+inv_gamma_lpdf(y[i] | shape,(shape-1)*scale[i]);
    } else
    {
    scale[i] = 0;
    log_lik[i] =0;
    count +=1;
    }
    if(bernoulli_rng(th[i])|| scale[i]<=0 || is_nan(scale[i])) yeval[i] =0; else yeval[i]=inv_gamma_rng(shape , scale[i]);
  }

  for (j in 1:NP) {

    if ((xpred[j] > tmin) && (xpred[j] < tmax)) th2[j]=inv_logit_shift(exp(-a+ nn*log(xpred[j]-tmin)+mm*log(tmax-xpred[j])),r,an); else th2[j]=1;
mu[j] = exp(-a+ nn*log(xpred[j]-tmin)+mm*log(tmax-xpred[j]));
if(bernoulli_rng(th2[j]) || mu[j]<=0 || is_nan(mu[j]))
{
  ypred[j]=0;
    mu[j] =0;
}else{
   ypred[j] =(inv_gamma_rng(shape,(shape-1)*mu[j]));
  }}
  alpha= exp(-a);
  dev = -2*sum(log_lik[]);
  tdmax = (nn*tmax+mm*tmin)/(nn+mm);
}"
    return(content)
  }
#' #' A string that contains the R-Stan code for Lactin model 
#' with Gaussian likelihood
#' @param funcode: no parameter is required
#' @return The output is the corresponding R-Stan code
#' @examples
#' output <- flactin_gauss();
#' @export
flactin_gauss=function(){
  content<-"functions { 
   vector flactin(vector y, vector theta,
 real[] x_r, int[] x_i){
 real x = y[1];
 real l = theta[1];
 real ro = theta[2];
 real a = theta[3];
 real del = theta[4];
 real val;
 vector[1] y_x;
 val = -l +  exp(ro*x) - a*exp(del*x);
 y_x[1] = val;
return y_x;                              
}
   real Sinv_lpdf(vector y, vector x, real shape, vector sc, real a, real l, real ro , real del) {
    real value;
    //real rat;
    int N;
    N=num_elements(y);
    value = shape * N * log(shape-1) - shape*N*l + shape*sum(exp(x*ro)) - shape*a*sum(exp(x*del)) - N*lgamma(shape) - (shape+1)*sum(log(y)) - sum(sc ./ y);
    return value;
  }
    real Sinvs_lpdf(real y, real sigma, real sc) {
    real value;
    real pi2=sqrt(2*pi());
    //real rat;
    //rat = shape * log(shape-1) + shape*N*l + shape*sum(exp(x*ro)) - shape*a*sum(exp(x*del));
    //value = shape *  log(shape-1) - shape*l + shape*exp(x*ro) - shape*a*exp(x*del) - lgamma(shape) - (shape+1)*(log(y)) - (sc / y);
    
    value = -log(pi2*sigma)-(1.0/2.0)*(1/sigma)*(y-sc)*(y-sc)*(1/sigma);
    return value;
  }
  }
data {
int<lower=0> N;
int<lower=0> NP;  
vector[NP] xpred;
vector[N] y;
vector[N] x;
vector[N] status;
}
transformed data{
vector[N] ey=exp(y);
real l1=min(x)-min(x)/2;
real l2=max(x)+max(x)/2;
vector[1] y_guess1;
vector[1] y_guess2;
real x_r[0];
int x_i[0];
int M=10;
//real st=1e-3;
y_guess1[1]= l1;
y_guess2[1]= l2;
}
parameters {
  //real <lower=0, upper=1> del; //lower 0
    real <lower=0, upper=1> del; //<lower=0>
  real <lower=0, upper=del>  ro; // free in previous case
  real <lower=0, upper=ro/del> a;
    //real <lower=0, upper=0.4> a;

  real <lower=0> l; //uper=1 previous value
real <lower=0>  sigmasq;
}
transformed parameters{
  real sigma;
 vector[4] theta;
 real  <lower=0> tmax;
 tmax = log(a)/(ro-del);
 theta[1]=l;
 theta[2]=ro;
 theta[3]=a;
 theta[4]=del;
  sigma = sqrt(sigmasq);

}
model {
vector[N] mu;
real count;
real cntr; // center of the temperature
real dst;
real ag;
real bg;
for (i in 1:N)
  {
    mu[i] = -l +  exp(ro*x[i]) - a*exp(del*x[i]); //exp((ro-del)*tmax) is equivalent to a
  y[i] ~ normal(mu[i], sigma);
  }
//target += gamma_lpdf(tmax|0.1, 0.01);
tmax ~ gamma(0.1, 0.01);
del ~ beta(1, 1);
l ~ gamma(0.1, 0.1);
ro ~ beta(1, 1);
sigmasq ~ inv_gamma(1e-3, 1e-3);

}
generated quantities { 
real log_lik[N];
//real log_lik2[N];
real yeval[N];
real sscale[N];
real ypred[NP];
real sscale2[NP];
real mu[NP];
real dev;
//real dev2;
real  <lower=0, upper=tmax> tdmax;
vector[1] xmin;
vector[1] xmax;
real s2cale=0;
real delta;
delta=1/del;
xmin = algebra_solver(flactin, y_guess1, theta, x_r, x_i, 1e-10, 1e+1, 1e+50);
xmax = algebra_solver(flactin, y_guess2, theta, x_r, x_i, 1e-10, 1e+1, 1e+50);

tdmax = tmax-log(ro/del)/(ro-del);

for (i in 1:N) {
  log_lik[i]=0;
   //log_lik2[i]=0;
     yeval[i]=0;

sscale[i] = -l +  exp(ro*x[i]) - a*exp(del*x[i]);
   
 if(!is_nan(sscale[i]) && !is_inf(sscale[i])){
   log_lik[i]=normal_lpdf( y[i] |sscale[i],sigma);
    yeval[i] = normal_rng(sscale[i],sigma);

}else 
{sscale[i]=0;
log_lik[i]=0;
yeval[i]=0;
}
  
}
for (j in 1:NP) {
sscale2[j] = -l +  exp(ro*xpred[j]) - a*exp(del*xpred[j]);  
if(is_nan(sscale2[j]))
{
  sscale2[j]=0;
  ypred[j]=0;
}else{
  ypred[j] = normal_rng(sscale2[j],sigma);
}
mu[j]= sscale2[j];
}
dev = -2*sum(log_lik[]);

} "
return(content)
}
#' #' A string that contains the R-Stan code for Lactin model 
#' with Inverse Gamma likelihood. No zeros considered 
#' @param funcode: no parameter is required
#' @return The output is the corresponding R-Stan code
#' @examples
#' output <- flactin_ig();
#' @export
flactin_ig=function(){
  content<-"functions { 
   vector flactin(vector y, vector theta,
 real[] x_r, int[] x_i){
 real x = y[1];
 real l = theta[1];
 real ro = theta[2];
 real a = theta[3];
 real del = theta[4];
 real val;
 vector[1] y_x;
 val = -l +  exp(ro*x) - a*exp(del*x);
 y_x[1] = val;
return y_x;                              
}
   real Sinv_lpdf(vector y, vector x, real shape, vector sc, real a, real l, real ro , real del) {
    real value;
    int N;
    N=num_elements(y);
    value = shape * N * log(shape-1) - shape*N*l + shape*sum(exp(x*ro)) - shape*a*sum(exp(x*del)) - N*lgamma(shape) - (shape+1)*sum(log(y)) - sum(sc ./ y);
    return value;
  }
    real Sinvs_lpdf(real y, real x, real shape, real sc, real a, real l, real ro , real del) {
    real value;
    value = shape *  log(shape-1) - shape*l + shape*exp(x*ro) - shape*a*exp(x*del) - lgamma(shape) - (shape+1)*(log(y)) - (sc / y);
    return value;
  }
    real Sin_lpdf(real y, real x, real shape, real sc, real a, real l, real ro , real del) {
    real value;
    value = shape *  log(shape-1) + shape*log(-l + exp(x*ro) - a*exp(x*del)) - lgamma(shape) - (shape+1)*(log(y)) - ((shape-1)*sc / y);
    return value;
  }  }
data {
  int<lower=0> N;
vector[N] y;
vector[N] x;
int<lower=0> NP;
vector[NP] xpred;

}
transformed data{
vector[N] ey=exp(y);
real l1=-min(x)-min(x)/2;
real l2=max(x)+max(x)/2;
vector[1] y_guess1;
vector[1] y_guess2;
real x_r[0];
int x_i[0];
y_guess1[1]= l1;
y_guess2[1]= l2;
}
parameters {

  real <lower=0, upper=1> del; //<lower=0>
  real <lower=0, upper=del>  ro; // free in previous case
   real  <lower=((log(ro)-log(del))/(ro-del))> tmax;
  real  <lower=-1, upper=1> l; 
  real <lower=2> shape;
}
transformed parameters{
 vector[4] theta;
 real <lower=0, upper=(ro/del)> a;
 a=exp((ro-del)*tmax);
 theta[1]=l;
 theta[2]=ro;
 theta[3]=a;
 theta[4]=del;
}
model {
vector[N] mu;
real count;

for (i in 1:N)
{
  mu[i] =  (-l+exp(ro*x[i])-a*exp(del*x[i]));

  if(is_nan(mu[i]) ) target += inv_gamma_lpdf( y[i] |shape, (shape-1)*mu[i]); else mu[i]=0;
 
 
}

//target += gamma_lpdf(tmax|0.1, 0.01);
tmax ~ gamma(0.1, 0.01);
del ~ beta(1, 1);
ro ~ beta(1, 1);
//l ~ normal(0, 10);
shape ~ gamma(.1, .01); 
}
generated quantities { 
real log_lik[N];
real yeval[N];
real mu[NP];
real ypred[NP];
real sscale2[NP];
real sscale[N];
real dev;
real  <lower=0, upper=tmax> tdmax;
vector[1] xmin;
vector[1] xmax;
real s2cale=0;
xmin = algebra_solver(flactin, y_guess1, theta, x_r, x_i, 1e-10, 1e+1, 1e+50);
xmax = algebra_solver(flactin, y_guess2, theta, x_r, x_i, 1e-10, 1e+1, 1e+50);
tdmax = tmax-log(ro/del)/(ro-del);

for (i in 1:N) {
sscale[i] = (-l+exp(ro*x[i])-a*exp(del*x[i]));

if (!is_inf(sscale[i]) && sscale[i]>0)

{        
  log_lik[i] = inv_gamma_lpdf( y[i] |shape, (shape-1)*sscale[i]);
  if (sscale[i]>0) yeval[i] = inv_gamma_rng(shape,(shape-1)*sscale[i]); else yeval[i]=0;


}else 
{
     yeval[i]=0;
     log_lik[i]=0;
     sscale[i]=0;
}
}
for (j in 1:NP) {
  sscale2[j] = (-l+exp(ro*xpred[j])-a*exp(del*xpred[j]));
 
if ( sscale2[j]>0 && (!is_inf(sscale2[j])) )
  {
    mu[j]=sscale2[j];
    ypred[j]=inv_gamma_rng(shape,(shape-1)*sscale2[j]);
  } else 
  {
    sscale2[j]=0;
    ypred[j]=0;
    mu[j]=0;
  }
}
dev = -2*sum(log_lik[]);
} "

return(content)
}
#' #' A string that contains the R-Stan code for Lactin model 
#' with zero inflated Inverse Gamma likelihood. Excess zeros considered 
#' @param funcode: no parameter is required
#' @return The output is the corresponding R-Stan code
#' @examples
#' output <- flactin_ziig();
#' @export
flactin_ziig=function(){
  content<-"functions {
  vector flactin(vector y, vector theta,
                 real[] x_r, int[] x_i){
    real x = y[1];
    real l = theta[1];
    real ro = theta[2];
    real a = theta[3];
    real del = theta[4];
    real val;
    vector[1] y_x;
    val = -l +  exp(ro*x) - a*exp(del*x);
    y_x[1] = val;
    return y_x;
  }
  real Sin_lpdf(real y, real x, real shape, real sc, real a, real l, real ro , real del) {
    real value;
    value = shape *  log(shape-1) + shape*log(-l + exp(x*ro) - a*exp(x*del)) - lgamma(shape) - (shape+1)*(log(y)) - ((shape-1)*sc / y);
    return value;
  }
  real inv_logit_shift(real x, real r, real an){
    real val;
    real test=an*(x-r);
      val=exp(- log_sum_exp(0,test));

    return val;
  }
}
data {
  int<lower=0> N;
  int<lower=0> NP;
  vector[NP] xpred;
  vector[N] y;
  vector[N] x;
  vector[N] status;
}
transformed data{
  real an=100;
  real r=5e-3; //r=3.006212e-04; // cause 0.05 is the median response and 0.1 is the max response.
  real l1=min(x)-min(x)/2;
  real l2=max(x)+max(x)/2;
  vector[1] y_guess1;
  vector[1] y_guess2;
  real x_r[0];
  int x_i[0];
  y_guess1[1]= l1;
  y_guess2[1]= l2;
}
parameters {
  real <lower=0, upper=1> del; //<lower=0>
    real <lower=0, upper=del>  ro; // free in previous case
    real  <lower=((log(ro)-log(del))/(ro-del))> tmax;
    real  l;
    real <lower=2> shape;

}
transformed parameters{
  real <lower=0, upper=(ro/del)> a;
  vector<lower=0, upper=1>[N] th;
  vector[4] theta;
  a=exp((ro-del)*tmax);
  theta[1]=l;
  theta[2]=ro;
  theta[3]=a;
  theta[4]=del;


  for (i in 1:N) {
    th[i]=inv_logit_shift((-l+exp(ro*x[i])-a*exp(del*x[i])),r,an);
  }
}
model {
  real log_other;
  real log_zero;
  vector[N]  ypred;
  for (i in 1:N)
  {

    if   (status[i]==0)   {
      ypred[i]=0;
      log_zero=log(th[i]);
      target += log_zero;
    } else{
      ypred[i] = (-l+exp(ro*x[i])-a*exp(del*x[i]));
      log_other=log(1-th[i]);
      target += log_other +Sin_lpdf(y[i]|x[i], shape, ypred[i], a, l,ro , del);

    }
  }

  tmax ~ gamma(0.1, 0.01);
  del ~ beta(1, 1);
  l ~ gamma(0.1, 0.1);
  ro ~ beta(1, 1);
  shape ~ gamma(0.1, 0.01);

}
generated quantities {
  real log_lik[N];
  real yeval[N];
  real scale[N];
  real dev;
  real tdmax;
  vector[1] xmin;
  vector[1] xmax;
  real ypred[NP];
  real sscale2[NP];
  real th2[NP];
  real mu[NP];

  real s2cale=0;
  xmin = algebra_solver(flactin, y_guess1, theta, x_r, x_i, 1e-10, 1e+1, 1e+10);
  xmax = algebra_solver(flactin, y_guess2, theta, x_r, x_i, 1e-10, 1e+1, 1e+10);

  for (i in 1:N) {
    scale[i] =  (-l+exp(ro*x[i])-a*exp(del*x[i]));
    if (status[i]==0)
    {
      if (!is_nan(th[i])) log_lik[i] = bernoulli_lpmf(1 | th[i]); else log_lik[i] = 0;
    }
    else if(!is_nan(scale[i]) && (!is_nan(th[i])) &&  scale[i]>0 && status[i]!=0)
    {
      log_lik[i] = bernoulli_lpmf(0 | th[i])+Sin_lpdf(y[i]|x[i], shape, scale[i], a, l,ro , del);
    } else
    {
      scale[i] = 0;
      log_lik[i] =0;
    }
    if(bernoulli_rng(th[i])|| scale[i]<=0 || is_nan(scale[i])) yeval[i] =0; else yeval[i]=inv_gamma_rng(shape , scale[i]);
  }
  for (j in 1:NP) {
    th2[j]=inv_logit_shift((-l+exp(ro*xpred[j])-a*exp(del*xpred[j])),r,an);
    sscale2[j] =(-l+exp(ro*xpred[j])-a*exp(del*xpred[j]));;

    if(bernoulli_rng(th2[j]) || sscale2[j]<=0 || is_nan(sscale2[j]))
    {
      ypred[j]=0;
      sscale2[j] =0;
      mu[j] =0;
    }else{
      ypred[j] =(inv_gamma_rng(shape,(shape-1)*sscale2[j]));
      mu[j] = sscale2[j];

    }}
  dev = -2*sum(log_lik[]);
  tdmax = tmax-log(ro/del)/(ro-del);
}"
  return(content)
}

#' An S4 class to represent a developmental rates object.
#'
#' @slot balance A length-one numeric vector
#' 
Account <- setClass("Account",
                    slots = list(balance = "numeric")
)
#' Fahrenheit conversion
#'
#' Convert degrees Fahrenheit temperatures to degrees Celsius
#' @param F_temp The temperature in degrees Fahrenheit
#' @return The temperature in degrees Celsius
#' @examples 
#' temp1 <- F_to_C(50);
#' temp2 <- F_to_C( c(50, 63, 23) );
#' @export
F_to_C <- function(F_temp){
  C_temp <- (F_temp - 32) * 5/9;
  return(C_temp);
}
#' Celsius conversion
#'
#' Convert degrees Celsius temperatures to degrees Fahrenheit
#' @param C_temp The temperature in degrees Celsius
#' @return The temperature in degrees Fahrenheit
#' @examples 
#' temp1 <- C_to_F(22);
#' temp2 <- C_to_F( c(-2, 12, 23) );
#' @export
C_to_F <- function(C_temp){
  F_temp <- (C_temp * 9/5) + 32;
  return(F_temp);
}
#' An S4 class to represent a developmental rates object.
#' Supresses output of a function
#' @param funcode: any function code
#' @return The output of the function without displaying it
#' @examples
#' output <- hush (print(1:100));
#' @export
hush=function(funcode){
  sink("NUL")
  tmp <- funcode
  sink()
  return(tmp)
}
#' Turn to uppercase the first letter of a string
#' @param x: any string
#' @return The string with first letter capital
#' @examples
#' output <- .simpleCap ("the quick red fox jumps over the lazy brown dog");
#' @export
.simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}
#' Create posterior predictions plot for Developmental Rate model
#' @param rdevel object that contains (est) s4 type Stanfit object, (summary) Stanfit summary of basic parameters, (diagnostics) Stanfit overal diagnostics report
#' @export
#' @seealso \code{\link[drin_hmc]{stan}}
#' @return Create posterior predictions plot for current model

drin_popp <- function(data, est, mtype="bieri",lik="gauss")
{
  yypred_mtype<-(as.matrix(est, pars="ypred"))
  ymeddata<-data.frame(x = data$xpred, y =   apply(yypred_mtype,2, median)) #for median
  yquant1data<-data.frame(x = data$xpred, y =   apply(yypred_mtype,2, function(x) quantile(x, c(0.025)) )) #for 1st qr
  yquant2data<-data.frame(x = data$xpred, y =   apply(yypred_mtype,2, function(x) quantile(x, c(0.975)) )) #for 2nd qr
  yquant3data<-data.frame(x = data$xpred, y =   apply(yypred_mtype,2, function(x) quantile(x, c(0.25)) )) #for 1st qr
  yquant4data<-data.frame(x = data$xpred, y =   apply(yypred_mtype,2, function(x) quantile(x, c(0.75)) )) #for 2nd qr
  
  
  #yspline_int <- as.data.frame(spline(ymeddata$x, ymeddata$y))
  
  ndata_mtype<-as.data.frame(data[2:4])
  yndata_mtype<-ymeddata
  # convert grouping variables to factor
  yndata_mtype$Legend<-as.factor(yndata_mtype$x)
  yndata_mtype$y1 <- yquant1data$y
  yndata_mtype$median <- ymeddata$y
  yndata_mtype$y2 <- yquant2data$y
  yndata_mtype$y3 <- yquant3data$y
  yndata_mtype$y4 <- yquant4data$y
  Legend<-factor(data$xpred) #the temperature levels
  colval<-colorRampPalette(c("royalblue","red"))
  #colval<-c("blue","#009933","#ffff00","#ff9933","#ff0000","#800000")
  ggplot2::scale_color_brewer(palette="Dark2")
  ggplot2::ggplot(data = yndata_mtype, aes(x = x, y = y)) + 
    #ylim(-0.75, 0.12)+
    ggplot2::geom_point(data=ndata_mtype, size = 3, aes(x = x, y = y, colour = "#000000"),alpha=1) +
    #    geom_errorbar(aes(ymin=y1, ymax=y2),width=.6,
    #                  position=position_dodge(0.5))+
    ggplot2::geom_line(aes(x = x, y = median, colour="median"), size=1)+
    ggplot2::geom_ribbon(aes(x = x,ymin=y1,ymax=y2),alpha=0.25, fill="#81DAF5")+
    ggplot2::geom_ribbon(aes(x = x,ymin=y3,ymax=y4),alpha=0.25, fill="#145fa7")+
    ggplot2::labs(
      x = expression("Temperature"~degree*C),
      y = expression("Developmental rate"),
      title = paste0("95% Cr.I. Posterior predictive plots based on",.simpleCap(mtype)," model") 
      #title = "Lactin model: Posterior Fitted Credible Intervals vs observed rates",
      #subtitle = "for Propylea quatuordecimpunctata"
    )+
    #subtitle = "for Tetranychus urticae mites")+
    #panel_bg(fill = "gray95", color = NA) +
    #grid_lines(color = "gray")+
    ggplot2::scale_colour_manual("Legend",labels = c("Data","Median"),
                        values = c("black","red"),
                        guide = guide_legend(override.aes = list(
                          linetype = c(rep("blank", 1),  "solid"),
                          shape = c(16, NA,NA),
                          fill=c("#81DAF5", "#81DAF5"))))+
    ggplot2::guides( color = guide_legend(
      order = 1,
      override.aes = list(
        color = c("black", "red"),
        fill  = c("white", "white"),
        linetype = c("blank", "solid"),
        shape = c(16, NA)))) +
    ggplot2::theme_bw() +
    # remove legend key border color & background
    ggplot2::theme(text = element_text(size=16),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.box.background = element_blank())+
    ggplot2::geom_point(aes(x = x, y = y, size = "95% Credible Interval", shape = NA)) +
    ggplot2::geom_point(aes(x = x, y = y, size = "50% Credible Interval", shape = NA)) +
    ggplot2::guides(size = guide_legend(NULL, 
                               order = 2,
                               override.aes = list(shape = c(15,15), 
                                                   color = c("#145fa7","#81DAF5"),
                                                   size = c(6,6))))+
    ggplot2::theme(legend.position = "none")
  
  
}
#' Developmental Rate using HMC
#'
#' Estimation of Ecological model parameters using HMC sampling
#' @param C_data data list including: y the response values, x the predictor values, status of censoring 0 if anthropod is not evolved ;
#' NP is the number predictions to be made. 
#' @param sparam vector indicating the number of chains, iterrations, burn_in and thinning in Rstan execution. Default is c(3,11000,10000,1);
#' @param mtype is the ecological model type between "bieri" (default), "briere", "analytis" and "lactin";
#' @param lik is the likelihood choice between "gauss" (default), "igamma" (inverse gamma) in rstan form;
#' @param prior is a list with parameter priors. If not declared default weekly informative;
#' @param init_list is the initial values of the parameters for the shake of the HMC algorithm. It should be a list or a function of a list 
#' @export
#' @seealso \code{\link[rstan]{stan}}
#' @return List that contains: (est) s4 type Stanfit object, (summary) Stanfit summary of basic parameters, (diagnostics) Stanfit overal diagnostics report
#' @examples \dontrun{
#' Do not run
#' DataBoot<-xls::read.xlsx(paste(bugname,".xlsx",sep=""), sheetIndex=1, header=FALSE, endRow = 1000)
#' DataBoot<-na.omit(DataBoot)
#' if (inum==1) {colnames(DataBoot)<- c("Temp","gr","status")} # if censored zero data exist
#' else {colnames(DataBoot)<- c("Temp","gr")} # if no zeros exist
#' y<-c(DataBoot$gr) #Observed days until development ocurrence
#' temp<-c(DataBoot$Temp) #create predictor for observed values (i.e. temperatures)
#' N <-dim(DataBoot)[1] #number of data
#' if (inum==1) status<-c(DataBoot$status) else status<-rep(1,N)
#' y<-1/y # create rensponse variable of observed developmental rates
#' if (inum==1) {data <- list(N = N, y = y, x=temp, status= status, t=10.0)} else {data <- list(N = N, y = y, x=temp, status= status, t=10.0)}  
#' data$NP<-round(2*(max(data$x)-min(data$x)),0) #number of predictions
#' data$xpred<-seq(min(data$x),max(data$x),length.out = data$NP) #generate NP predictor values for model predictions
#' modbieri_stan<-drin_hmc(data,c(1,11000,10000,1)) #run the developmental 
#' modbieri_stan<-drin_hmc(data,c(1,11000,10000,1),"bieri","gauss")
#' }
drin_hmc <- function(c_data, sparam=c(4,11000,10000,1), mtype="bieri",lik="gauss",prior="NULL",init_list="NULL",predplot=FALSE,...)
  {
  if (is.null(c_data$N)) c_data$N<-length(c_data$x)
  #if (is.null(c_data$status)) c_data$status<-(rep(1,length(c_data$x)))
  c_data$status<-(rep(1,length(c_data$x)))
  if (is.null(c_data$t)) c_data$t<-10
  if (is.null(c_data$NP) && is.null(c_data$xpred)) c_data$NP<-round(2*(max(c_data$x)-min(c_data$x)),0) else if (!is.null(c_data$xpred)) c_data$NP<-length(c_data$xpred) 
  if (is.null(c_data$xpred)) c_data$xpred<-seq(min(c_data$x),max(c_data$x),length.out = c_data$NP)
  if (length(which(c_data$y==0))>0) zeropresent<- TRUE else zeropresent<- FALSE
  if (zeropresent) c_data$status[which(c_data$y==0)]<-0 else c_data$status<-(rep(1,length(c_data$x)))
  
  if (!is.null(sparam[1])) nch<-sparam[1] else nch<-4
  if (!is.null(sparam[2])) iter<-sparam[2] else iter<-11000
  if (!is.null(sparam[3])) bin<-sparam[3] else bin<-10000
  if (!is.null(sparam[4])) nth<-sparam[4] else nth<-1
  
  inittmin<-rgamma(1,shape=100,rate=10)

if (lik=="igamma") {
  switch(mtype,
         bieri = {if(zeropresent) cfile <- fbieri_ziig() else cfile <- fbieri_ig()
         brfun_bieri_ig<- function (num_id=1) {list( a=rgamma(1,shape=1e+5,rate=1e+6),
                                                     tmin=inittmin-rbeta(1,1,1),
                                                     b=1+runif(1,0.001,0.1),
                                                     tmax=inittmin+rgamma(1,shape=200,rate=10),
                                                     shape=runif(1,300,500))}
         brfun<-brfun_bieri_ig},
         briere = {if(zeropresent) cfile <- fbriere_ziig() else cfile <- fbriere_ig()
         brfun_briere_ig<- function (num_id=1)  {list(
           a=rgamma(1,shape=1e+2,rate=1e+1),
           tmin=inittmin-rbeta(1,1,1),
           tmax=inittmin+rgamma(1,shape=200,rate=10), 
           shape=100+10*runif(1,num_id,10))}
         brfun<-brfun_briere_ig},
         analytis = {if(zeropresent) cfile <- fanalytis_ziig() else cfile <- fanalytis_ig()
         brfun_analytis_ig<- function (num_id=1) {list(a=rgamma(1,shape=1e+2,rate=1e+1), 
                                                       tmin=max(inittmin-rbeta(1,1,1),4),
                                                       tmax=max(inittmin,4)+rgamma(1,shape=200,rate=10), 
                                                       nn=1+rbeta(1,shape1=1,shape2=1), 
                                                       mm=0.07*rbeta(1,shape1=1,shape2=1), 
                                                       shape=100*runif(1,num_id,10))}
         brfun<-brfun_analytis_ig},
         lactin = {if(zeropresent) cfile <- flactin_ziig() else cfile <- flactin_ig()
         del=runif(1,0,0.2)
         th=rbeta(1,1,1)
         ro=del-(1/100)*rbeta(1,shape1=1,shape2=1)
         brfun_lactin_ig<- function (num_id=1) {list(tmax=38+rbeta(1,1,1),
                                                     a=((ro/del))*rbeta(1,1,1),
                                                     l=-1+2*rbeta(1,1,1), 
                                                     del=del, 
                                                     ro=ro, 
                                                     shape=200+rgamma(1,0.1,0.01)
                                                     #th=rbeta(1,1,1),
                                                     #r=rgamma(1,1e-3,1e-1)
         )}
         brfun<-brfun_lactin_ig})
#  gsub("gauss", "ig", cfile
  #gsub("g", "ig", init)#substitude gauss to inverse gamma if selected
} else{
  switch(mtype,
         bieri = {cfile <- fbieri_gauss()
         ttmax=rgamma(1,50,1)
         brfun_bieri_g<- function (num_id=1) {list(st=rgamma(1,10,1),
                                                   a=rbeta(1,1,1)/100,
                                                   tmin=ttmax*rbeta(1,1,1),
                                                   b=1+runif(1,0.5,1), tmax=ttmax,
                                                   sigmasq=rgamma(1,shape=10,rate=10))}
         brfun<-brfun_bieri_g},
         briere = {cfile <- fbriere_gauss()
         brfun_briere_g<- function (num_id=1) {list(st=rgamma(1,1,1e-1),
                                          a=runif(1,5,15), 
                                          tmin=inittmin-rbeta(1,1,1),
                                          tmax=inittmin+rgamma(1,shape=200,rate=10), 
                                          sigmasq=runif(1,0,num_id))}
         brfun<-brfun_briere_g},
         analytis = {cfile <- fanalytis_gauss()
         brfun_analytis_g<- function (num_id=1) {list(st=rgamma(1,2,1),
                                                      a=rgamma(1,shape=1e+2,rate=1e+1), 
                                                      tmin=inittmin-rbeta(1,1,1),
                                                      tmax=inittmin+rgamma(1,shape=200,rate=10),                                                       nn=1+rbeta(1,shape1=1,shape2=1), 
                                                      mm=0.5*rbeta(1,shape1=1,shape2=1), 
                                                      sigmasq=runif(1,num_id,10))}
         brfun<-brfun_analytis_g},
         lactin = {cfile <- flactin_gauss()
         del=runif(1,0.5,1)
         ro=del*rbeta(1,shape1=1,shape2=1)
         brfun_lactin_g<- function (num_id=1) {list(st=rgamma(1,2,1),
                                                    a=((ro/del))*rbeta(1,1,1),
                                                    l=1+rgamma(1,1,2), 
                                                    del=del, 
                                                    ro=ro, 
                                                    sigmasq=runif(1,1,num_id*10))}
      
         brfun<-brfun_lactin_g})}

if ((!is.null(prior)) && (!is.null(names(prior)))) 
  for (i in 1:length(prior)) 
{
  nm<-names(prior)[i]
  nm_value<-paste0(nm," ~ ",prior[i],";")
  nm_pattern<-paste0("(",nm," ~)(.+)?[;]") # search name of parameter
  
  #nm_pattern<-"(del ~)(.+)?[;]"
  #nm_value<-paste0("del"," ~ ","normal(0,100)",";")
  #str(unlist(str_extract_all(cfile, nm_pattern)))
  
  cfile<-gsub(unlist(stringr::str_extract_all(cfile, nm_pattern)), nm_value, cfile, fixed = TRUE)
  #unlist(str_extract_all(cnew, nm_pattern))
  
  }
if ((!is.null(init_list)) && (!is.null(names(init_list)))) brfun<-init_list
  
  
  #for (j in 1:length(init_list)) 
  #{
    #rnm<-names(init_list)[j]
    #rnm_value<-paste0(rnm,"=",init_list[j],";")
    #rnm_pattern<-paste0("(",rnm,"=)(.+)?[;]") # search name of parameter
    
    #rnm<-"del"
    #rnm_value<-expression(rbeta(1,1,1))
    
    
    #brlist<-as.list(body(brfun)[[2]])
    #brlist[[rnm]]<-rnm_value
    #brtext<-paste(names(brlist),brlist,sep="=",collapse="," )
    #brtext<-sub(".*list,","",brtext)
    #sub(".*=", "", brtext)
    
    
    #eqstr<-strsplit(brtext, "=")
    
    
    #str_extract_all(brtext, "[a-z]+(.=)")
    
    #str_extract_all(brtext, ".+[=]")
    
    
    
   # sub("\\=.*", "", brtext)
    
  #  gsub(".*:", "", x)
  #  nmlist<-gsub(".*=","",brtext)
    
    
    
   # brfun<-function (num_id=1) {}
    
  #  as.list(strsplit(brtext, "),")[[1]])
    
  #  eval(parse(text=paste('list(', brtext, ')')))
  #  as.list(scan(text = brtext, what = "", sep = ","))
  #  as.list(el(strsplit(brtext, "),")))
    
   # stringr::str_split(string = brtext, pattern = "[)],")
    
  
  #    is.vector(brlist)
    
   # brlist[[2]][[3]]
  #body(brfun)[[2]] <- substitute(del <- 2)
    
    
  #  rnm_pattern<-"(del=)(.+)?[,]"
  #  rnm_pattern<-"(del=)(.+)"
  #  rnm_value<-paste0("del","=","normal(1,0,100)",",")
  #  str_extract_all(brlist[[2]], rnm_pattern)

    
   # rnm_check<-paste0(rnm,"=")
    #nm_check<-paste0("del"," ~ ")
  #  stringr::str_replace(brfun, unlist(stringr::str_extract_all(brfun, rnm_check)), rnm_value)
  #}


inits<-lapply(1:nch, function(id) brfun(num_id=id))
est <- rstan::stan(model_code = cfile, init = inits, save_dso=FALSE, control = list(max_treedepth = 11,adapt_delta = 0.99), algorithm = "NUTS", data = c_data, iter = iter, chains = nch, warmup=bin, thin=nth,verbose = TRUE,...)

output <- list(est=est)
  
parnames<- names(rstan::extract(est))
param<-parnames[which(parnames %in% c("a","b","mm","nn","ro","del","xmin","tmin","tdmax","tmax","xmax","shape","sigma","dev","lp__"))]

##est$call <- match.call()
#output$summary <- as.data.frame(rstan::summary(est,pars=param)$summary)
#output$diagnostics <- capture.output(rstan::check_hmc_diagnostics(est), type = "message", split = FALSE)

output$summary <- as.data.frame(rstan::summary(est,pars=param)$summary)
output$diagnostics <- hush(capture.output(rstan::check_hmc_diagnostics(est), type = "message", split = FALSE))
output$data <- c_data
output$mtype <- mtype
output$lik <- lik
print(output$summary)
if (predplot==TRUE) drin_popp(c_data, est, mtype=mtype,lik=lik) #plot posterior predictions
##class(est) <- "drin_hmc"
##est
return(output)
}





