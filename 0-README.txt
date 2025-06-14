The supplementary materials contain the matlab code and dataset in paper: "On Modeling of Multiplicative Measurement Errors for Multivariate Degradation Data".

The details of the supplementary materials are as follows:


(1) main_bearing.m

- Description: 
- This file is used to analyze the degradation data from bearing 1 and bearing 2.
- See details of the analysis in supplementary in "Bearings on high-speed trains".

- Input:
- bearing_1.mat and bearing_2.mat: degradation data on high-speed train from Bearing 1 and Bearing 2.

- Output:
- parameter_estimation_proposed_model_bearing.mat: Parameter estimates at each iteration of the MCEM algorithm.


(2) main_mcem_simulation.m

- Description:
- The file conducts one replication of MCEM algorithm in estimating model parameters in the proposed Wiener process.
- See details of the proposed multivariate Wiener process model in Equation (4) in the paper.
- See details of the MCEM algorithm in Section "Wiener processes under repeated tests" in the paper.
- The measurement scheme in the degradation test is repeated measurement.

- Functions required:
- Wiener_degradation_simulation_3d: simulation degradation data.
- Wiener_MCEM_initiation_3d: find initial values in model parameters.
- Wiener_MCEM_3d: MCEM algorithm.


(3) Wiener_degradation_simulation_3d.m 

- Description:
- This function simulates the degradation data using the proposed Wiener process model.
- See details of the proposed multivariate Wiener process model in Equation (4) in the paper.
- The measurement scheme in the degradation test is repeated measurement.

- Arguments:
- degra_params: the setting of degradation data.
- - degra_params.n: batch size.
- - degra_params.K: sample size in each batch.
- - degra_params.m: number of measurements for each unit.
- - degra_params.d: dimension of multivariate degradation.
- - degra_params.delta_t: time interval of measurement times.
- - degra_params.v: a vector, degradation rate in multivariate Wiener process.
- - degra_params.sigma2: a matrix, covariance matrix in multivairate Wiener process.
- - degra_params.eta: both the shape and rate parameters in random error.
- debug_flag: the setting of whether show the simulated data.

- Output:
- Y: degradation measurements.
- t: measurement times.
- dY: degradation increments.
- dt: time intervals.
- zeta: measurement errors.


(4) Wiener_MCEM_initiation_3d.m

- Description:
- This function initializes the parameters in the proposed Wiener process.
- See details of the proposed multivairate Wiener process model in Equation (4) in the paper.

- Arguments:
- dY: increments in degradation measurements.
- dt: increments in measurement times.
- n: number of batches.
- K: number of units in each batch.
- m: number of measurements for each unit.

- Output:
- v1_ini: initilization of drift parameter v1 for PC1.
- v2_ini: initilization of drift parameter v2 for PC2.
- v3_ini: initilization of drift parameter v3 for PC3.
- v_ini: vector of [v1_ini,v2_ini,v3_ini]'.
- sigma2_ini: initilization of the covariance matrix in multivairate Wiener process.
- eta_ini: initilization of the shape parameter eta in gamma distribution.
- zeta_ini: initilization of the measurement errors.


(5) Wiener_MCEM_3d.m

- Description:
- this function estimates the model parameters in Wiener process with 3 PCs using MCEM algorithm.
- See details of the proposed multivariate Wiener process model in Equation (4) in the paper.
- See details of the MCEM algorithm in Section "Wiener processes under repeated tests" in the paper.
- The measurement scheme in the degradation test is repeated measurement.

- Usage:
- [opts] = Wiener_MCEM_3d(dt,dY,n,m,d,MCEM_set_params,K,v_ini,sigma2_ini,eta_ini,zeta_ini)

- Arguments:
- dt: time intervals.
- dY: degradation increments.
- n: batch size.
- K: sample size in each batch.
- m: number of measurements for each unit.
- d: dimensions of multivariate degradation.
- v_ini: initiation for parameter v with shape d*1.
- sigma2_ini: initiation for parameter Sigma.
- eta_ini: initiation for parameter eta.
- zeta_ini: initiation for measurement errors.
- MCEM_set_params: settings in MCEM algorithm.
- - MCEM_set_params.B: number of MCMC samples in Algorithm 1.
- - MCEM_set_params.L: number of iterations in MCEM algorithm.
- - MCEM_set_params.tol: tolerance in MCEM algorithm.
- - MCEM_set_params.mcmc_sigma_c: tuning parameter in MCMC in Algorithm 1.
- - MCEM_set_params.mcmc_burn_percent: early B_0 samples are removed as burn-in.
- - MCEM_set_params.mcmc_interval: select one sample in every u samples in Algorithm 1.

- Outputs:
- opts: a struct saving all outputs.
- - opts.eta_est: estimates of eta.
- - opts.v_est: estimates of v, i.e., the degradation rates.
- - opts.v1_est: estimates of v1, i.e., degradation rate of the first PC.
- - opts.v2_est: estimates of v2, i.e., degradation rate of the second PC.
- - opts.v3_est: estimates of v3, i.e., degradation rate of the third PC.
- - opts.sigma2_est: estimates of sigma.
- - opts.sigma2_11_est: estimates of sigma_11.
- - opts.sigma2_12_est: estimates of sigma_12.
- - opts.sigma2_13_est: estimates of sigma_13.
- - opts.sigma2_22_est: estimates of sigma_22.
- - opts.sigma2_23_est: estimates of sigma_23.
- - opts.sigma2_33_est: estimates of sigma_33.
- - opts.zeta_est: estimates of varsigma, i.e., the measurement error.



