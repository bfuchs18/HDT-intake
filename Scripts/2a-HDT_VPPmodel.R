######################################################
###### FULL IGT_VPP CODE FROM AHN'S GROUP ######
######################################################

graphics.off()
rm(list=ls(all=TRUE))

install.packages("hBayesDM")
library(hBayesDM)
library(rstan)
library(StanHeaders)
library(loo)

# library(rstudioapi)
# current_path <- getSourceEditorContext()$path

# setwd(dirname(current_path))

# To see how long computations take
startTime <- Sys.time()


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

data = "HDT_raw-task-data.txt"
inits = "fixed"
indPars = "mean"

rawdata <- read.table(data, header = T, sep = "\t")

# Remove rows containing NAs 
NA_rows_all = which(is.na(rawdata), arr.ind = T)  # rows with NAs
NA_rows = unique(NA_rows_all[, "row"])
if (length(NA_rows) > 0) {
  rawdata = rawdata[-NA_rows,]
  cat("The number of rows with NAs = ", length(NA_rows), ". They are removed prior to modeling the data. \n", sep = "")
}


# Individual Subjects
subjList <- unique(rawdata[,"subjID"])  # list of subjects x blocks
numSubjs <- length(subjList)  # number of subjects

# Specify the number of parameters and parameters of interest
numPars <- 8
POI     <- c("mu_A", "mu_alpha", "mu_cons", "mu_lambda", "mu_epP", "mu_epN", "mu_K", "mu_w",
             "sigma",
             "A", "alpha", "cons", "lambda", "epP", "epN", "K", "w",
             "log_lik", "y_pred")

modelName <- "igt_vpp"

################################################################################
# THE DATA.  ###################################################################
################################################################################

Tsubj <- as.vector(rep(0, numSubjs)) # number of trials for each subject

for (sIdx in 1:numSubjs)  {
  curSubj     <- subjList[sIdx]
  Tsubj[sIdx] <- sum(rawdata$subjID == curSubj)  # Tsubj[N]
}

maxTrials <- max(Tsubj)

# Information for user continued
cat(" # of (max) trials per subject = ", maxTrials, "\n\n")

RLmatrix <- array(0, c(numSubjs, maxTrials))
Ydata    <- array(-1, c(numSubjs, maxTrials))

#RLmatrix gets filled in with the NetScore for each trial??
#Ydata gets filled in withe the door choice for each trial

for (subjIdx in 1:numSubjs)   {
  #number of trials for each subj.
  useTrials                      <- Tsubj[subjIdx]
  currID                         <- subjList[subjIdx]
  rawdata_curSubj                <- subset(rawdata, rawdata$subjID == currID)
  RLmatrix[subjIdx, 1:useTrials] <- rawdata_curSubj[, "gain"] - 1 * abs(rawdata_curSubj[, "loss"])
  
  for (tIdx in 1:useTrials) {
    Y_t                     <- rawdata_curSubj[tIdx, "choice"] # chosen Y on trial "t"
    Ydata[subjIdx , tIdx] <- Y_t
  }
}

payscale = .25

dataList <- list(
  N       = numSubjs,
  T       = maxTrials,
  Tsubj   = Tsubj,
  outcome = RLmatrix * payscale,
  choice  = Ydata
)

# inits
if (inits[1] != "random") {
  if (inits[1] == "fixed") {
    inits_fixed = c(0.5, 0.5, 1.0, 1.0, 0, 0, 0.5, 0.5)
  } else {
    if (length(inits) == numPars) {
      inits_fixed = inits
    } else {
      stop("Check your inital values!")
    }
  }
  genInitList <- function() {
    list(
      mu_p      = c(qnorm(inits_fixed[1]), qnorm(inits_fixed[2]), qnorm(inits_fixed[3] / 5), qnorm(inits_fixed[4] / 10),
                    inits_fixed[5], inits_fixed[6], qnorm(inits_fixed[7]), qnorm(inits_fixed[8])),
      sigma     = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
      A_pr      = rep(qnorm(inits_fixed[1]), numSubjs),
      alpha_pr  = rep(qnorm(inits_fixed[2]), numSubjs),
      cons_pr   = rep(qnorm(inits_fixed[3]/5), numSubjs),
      lambda_pr = rep(qnorm(inits_fixed[4]/10), numSubjs),
      epP_pr    = rep(inits_fixed[5], numSubjs),
      epN_pr    = rep(inits_fixed[6], numSubjs),
      K_pr      = rep(qnorm(inits_fixed[7]), numSubjs),
      w_pr      = rep(qnorm(inits_fixed[8]), numSubjs)
    )
  }
} else {
  genInitList <- "random"
}


###################################################
######## CODE FOR STAN MODEL (FROM GITHUb) ########
###################################################

modelString = "
data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsubj[N];
  real outcome[N, T];
  int choice[N, T];
}

transformed data {
  vector[4] initV;
  initV  = rep_vector(0.0, 4);
}

parameters {
  // Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters
  vector[8] mu_p;
  vector<lower=0>[8] sigma;
  
  // Subject-level raw parameters (for Matt trick)
  vector[N] A_pr;
  vector[N] alpha_pr;
  vector[N] cons_pr;
  vector[N] lambda_pr;
  vector[N] epP_pr;
  vector[N] epN_pr;
  vector[N] K_pr;
  vector[N] w_pr;
}

transformed parameters {
  // Transform subject-level raw parameters
  vector<lower=0, upper=1>[N]  A;
  vector<lower=0, upper=2>[N]  alpha;
  vector<lower=0, upper=5>[N]  cons;
  vector<lower=0, upper=10>[N] lambda;
  vector[N] epP;
  vector[N] epN;
  vector<lower=0, upper=1>[N] K;
  vector<lower=0, upper=1>[N] w;
  
  for (i in 1:N) {
    A[i]      = Phi_approx(mu_p[1] + sigma[1] * A_pr[i]);
    alpha[i]  = Phi_approx(mu_p[2] + sigma[2] * alpha_pr[i]) * 2;
    cons[i]   = Phi_approx(mu_p[3] + sigma[3] * cons_pr[i]) * 5;
    lambda[i] = Phi_approx(mu_p[4] + sigma[4] * lambda_pr[i]) * 10;
    K[i]      = Phi_approx(mu_p[7] + sigma[7] * K_pr[i]);
    w[i]      = Phi_approx(mu_p[8] + sigma[8] * w_pr[i]);
  }
  epP = mu_p[5] + sigma[5] * epP_pr;
  epN = mu_p[6] + sigma[6] * epN_pr;
}

model {
  // Hyperparameters
  mu_p[1]  ~ normal(0, 1.0);
  mu_p[2]  ~ normal(0, 1.0);
  mu_p[3]  ~ normal(0, 1.0);
  mu_p[4]  ~ normal(0, 1.0);
  mu_p[5]  ~ normal(0, 10.0);
  mu_p[6]  ~ normal(0, 10.0);
  mu_p[7]  ~ normal(0, 1.0);
  mu_p[8]  ~ normal(0, 1.0);
  sigma ~ cauchy(0, 5);
  
  // individual parameters
  A_pr      ~ normal(0, 1.0);
  alpha_pr  ~ normal(0, 1.0);
  cons_pr   ~ normal(0, 1.0);
  lambda_pr ~ normal(0, 1.0);
  epP_pr    ~ normal(0, 1.0);
  epN_pr    ~ normal(0, 1.0);
  K_pr      ~ normal(0, 1.0);
  w_pr      ~ normal(0, 1.0);
  
  for (i in 1:N) {
    // Define values
    vector[4] ev;
    vector[4] p_next;
    vector[4] str;
    vector[4] pers;   // perseverance
    vector[4] V;   // weighted sum of ev and pers
    
    real curUtil;     // utility of curFb
    real theta;       // theta = 3^c - 1
    
    // Initialize values
    theta = pow(3, cons[i]) -1;
    ev    = initV; // initial ev values
    pers  = initV; // initial pers values
    V     = initV;
    
    for (t in 1:Tsubj[i]) {
      // softmax choice
      choice[i, t] ~ categorical_logit(theta * V);
      
      // perseverance decay
      pers = pers * K[i]; // decay
      
      if (outcome[i, t] >= 0) {  // x(t) >= 0
        curUtil = pow(outcome[i, t], alpha[i]);
        pers[choice[i, t]] = pers[choice[i, t]] + epP[i];  // perseverance term
      } else {                  // x(t) < 0
        curUtil = -1 * lambda[i] * pow(-1 * outcome[i, t], alpha[i]);
        pers[choice[i, t]] = pers[choice[i, t]] + epN[i];  // perseverance term
      }
      
      ev[choice[i, t]] = ev[choice[i, t]] + A[i] * (curUtil - ev[choice[i, t]]);
      // calculate V
      V = w[i] * ev + (1-w[i]) * pers;
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1>  mu_A;
  real<lower=0, upper=2>  mu_alpha;
  real<lower=0, upper=5>  mu_cons;
  real<lower=0, upper=10> mu_lambda;
  real mu_epP;
  real mu_epN;
  real<lower=0, upper=1> mu_K;
  real<lower=0, upper=1> mu_w;
  
  // For log likelihood calculation
  real log_lik[N, T];
  
  // For posterior predictive check
  real y_pred[N, T];
  
  // Set all posterior predictions to 0 (avoids NULL values)
  for (i in 1:N) {
    for (t in 1:T) {
      y_pred[i, t] = -1;
    }
  }
  
  mu_A      = Phi_approx(mu_p[1]);
  mu_alpha  = Phi_approx(mu_p[2]) * 2;
  mu_cons   = Phi_approx(mu_p[3]) * 5;
  mu_lambda = Phi_approx(mu_p[4]) * 10;
  mu_epP    = mu_p[5];
  mu_epN    = mu_p[6];
  mu_K      = Phi_approx(mu_p[7]);
  mu_w      = Phi_approx(mu_p[8]);
  
  { // local section, this saves time and space
    for (i in 1:N) {
      // Define values
      vector[4] ev;
      vector[4] p_next;
      vector[4] str;
      vector[4] pers;   // perseverance
      vector[4] V;   // weighted sum of ev and pers
      
      real curUtil;     // utility of curFb
      real theta;       // theta = 3^c - 1
      
      // Initialize values
      theta      = pow(3, cons[i]) -1;
      ev         = initV; // initial ev values
      pers       = initV; // initial pers values
      V          = initV;
      
      for (t in 1:Tsubj[i]) {
        // softmax choice
        log_lik[i, t] = categorical_logit_lpmf(choice[i, t] | theta * V);
        
        // generate posterior prediction for current trial
        y_pred[i, t] = categorical_rng(softmax(theta * V));
        
        // perseverance decay
        pers = pers * K[i]; // decay
        
        if (outcome[i, t] >= 0) {  // x(t) >= 0
          curUtil = pow(outcome[i, t], alpha[i]);
          pers[choice[i, t]] = pers[choice[i, t]] + epP[i];  // perseverance term
        } else {                  // x(t) < 0
          curUtil = -1 * lambda[i] * pow(-1 * outcome[i, t], alpha[i]);
          pers[choice[i, t]] = pers[choice[i, t]] + epN[i];  // perseverance term
        }
        
        ev[choice[i, t]] = ev[choice[i, t]] + A[i] * (curUtil - ev[choice[i, t]]);
        // calculate V
        V = w[i] * ev + (1-w[i]) * pers;
      }
    }
  }
}
"

writeLines(modelString , con = "HDTkids_vpp.stan") 

fit = stan(file   = "HDTkids_vpp.stan",
           data   = dataList,
           pars   = POI,
           warmup = 1000,
           init   = genInitList,
           iter   = 5000,
           chains = 2,
           thin   = 1,
           control = list(adapt_delta   = 0.99,
                          max_treedepth = 25,
                          stepsize      = 1))

LLVPP = extract_log_lik(fit, parameter_name = "log_lik")


looVPP  <- loo(LLVPP)

pareto_k_table(looVPP)

save.image(sprintf("VPPKids%s.Rdata", Sys.Date()))
