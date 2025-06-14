%{

- Description:
- The file conducts one replication of MCEM algorithm in estimating model parameters in Wiener process model.
- See details of the proposed multivariate Wiener process model in Equation (4) in the paper.
- See details of the MCEM algorithm in Section "Wiener processes under repeated tests" in the paper.
- The measurement scheme in the degradation test is repeated measurement.

- Functions required:
- Wiener_degradation_simulation_3d: simulation degradation data.
- Wiener_MCEM_initiation_3d: find initial values in model parameters.
- Wiener_MCEM_3d: MCEM algorithm.

- Notation:
- in the code       | in the paper    | Explanations
- n                 | n               | number of batches in degradation test
- K                 | K               | number of units in each batch i
- m                 | m               | the total number of measurements in batch i
- d                 | d               | number of PCs in each unit
- delta_t           | delta_t         | time interval between measurements
- v                 | v               | degradation rates in Wiener process
- sigma2            | Sigma           | covariance matrix in Wiener process
- eta               | eta             | both the shape and rate parameters in gamma distribution
- _est              | \               | estimation of the parameters
- Y                 | Y               | degradation levels of the d PCs
- t                 | t               | measurement times
- B                 | B_0+Mu          | number of MCMC samples  
- L                 | \               | number of iterations in MCEM algorithm
- tol               | \               | tolerance in MCEM algorithm
- mcmc_sigma_c      | sigma_c         | tuning parameter in MCMC 
- mcmc_burn_percent | B_0/B           | early B_0 samples are removed as burn-in
- mcmc_interval     | u               | select one sample in every u samples  

%}

clc;
clear all;
close all;

% batch size
n = 7;
% sample size in each batch
K = 3;
% number of measurements
m = 20;
% number of degradation characteristics
d = 3;
delta_t = 40/m;
% parameters setting
v = [2,4,6]';
sigma2_11 = 1.7;
sigma2_22 = 3.2;
sigma2_33 = 6.3;
rho_12 = 0.63;
rho_13 = 0.61;
rho_23 = 0.65;
sigma2_12 = rho_12*sqrt(sigma2_11*sigma2_22);
sigma2_13 = rho_13*sqrt(sigma2_11*sigma2_33);
sigma2_21 = sigma2_12;
sigma2_23 = rho_23*sqrt(sigma2_22*sigma2_33);
sigma2_31 = sigma2_13;
sigma2_32 = sigma2_23;
% set covariance matrix
sigma2 = [sigma2_11 sigma2_12 sigma2_13;
    sigma2_21 sigma2_22 sigma2_23;
    sigma2_31 sigma2_32 sigma2_33];
eta = 50;
% save the pre_set parameters in a struct
degra_params = struct();
degra_params.n = n;
degra_params.K = K;
degra_params.m = m;
degra_params.d = d;
degra_params.delta_t = delta_t;
degra_params.v = v;
degra_params.sigma2 = sigma2;
degra_params.rho_12 = rho_12;
degra_params.rho_13 = rho_13;
degra_params.rho_23 = rho_23;
degra_params.eta = eta;

% degradation data simulation
debug_flag = 0;
[Y,t,dY,dt,zeta] = Wiener_degradation_simulation_3d(degra_params,debug_flag);

% Initialization for MCEM
[v1_ini,v2_ini,v3_ini,v_ini,sigma2_ini,eta_ini,zeta_ini] = Wiener_MCEM_initiation_3d(dY,dt,n,K,m);

% MCEM algorithm
show_bias = 1;
B = 1000;
L = 50;
tol = 0.05;
mcmc_sigma_c = 2;
mcmc_burn_percent = 0.05;
mcmc_interval = 20;
MCEM_set_params.B = B;
MCEM_set_params.L = L;
MCEM_set_params.tol = tol;
MCEM_set_params.mcmc_sigma_c = mcmc_sigma_c;
MCEM_set_params.mcmc_burn_percent = mcmc_burn_percent;
MCEM_set_params.mcmc_interval = mcmc_interval;
[MCEM_opts] = Wiener_MCEM_3d(dt,dY,n,m,d,MCEM_set_params,K,v_ini,sigma2_ini,eta_ini,zeta_ini);
% save the MCEM results
eta_est = MCEM_opts.eta_est;
v_est = MCEM_opts.v_est;
v1_est = MCEM_opts.v1_est;
v2_est = MCEM_opts.v2_est;
v3_est = MCEM_opts.v3_est;
sigma2_est = MCEM_opts.sigma2_est;
sigma2_11_est = MCEM_opts.sigma2_11_est;
sigma2_12_est = MCEM_opts.sigma2_12_est;
sigma2_13_est = MCEM_opts.sigma2_13_est;
sigma2_22_est = MCEM_opts.sigma2_22_est;
sigma2_23_est  =MCEM_opts.sigma2_23_est;
sigma2_33_est = MCEM_opts.sigma2_33_est;
zeta_est = MCEM_opts.zeta_est;
zeta_now = MCEM_opts.zeta_now;