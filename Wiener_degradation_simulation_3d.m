function [Y,t,dY,dt,zeta] = Wiener_degradation_simulation_3d(degra_params,debug_flag)
%{

- Description:
- This function simulates the degradation data Using proposed Wiener process model.
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


- Model:
- The baseline degradation process is Wiener process: 
- X_{ik}(t) = vt+sigma2^(1/2)*B(t)
- The degradation observations are:
- Y_{ijk} = varsigma_{ij}*X_{ik}(t_{ij})
- varsigma_{ij} follows gamma distribution gamma(eta,eta)
- mean of varsigma : 1
- variance of varsigma : 1/eta
- we use "zeta" to denote "varsigma" in the code, for convenience

%}

n = degra_params.n;
K = degra_params.K;
m = degra_params.m;
d = degra_params.d;
delta_t = degra_params.delta_t;
v = degra_params.v;
sigma2 = degra_params.sigma2;
eta = degra_params.eta;
% random seed
rng(26,'philox');
seed = rng;
% set measurement times
for i = 1:1:n
    dt(i,1,1:K) = delta_t;
    t(i,1,1:K) = delta_t;
    for j = 2:1:m
        for k = 1:1:K
            t(i,j,k) = delta_t*j;
            dt(i,j,k) = t(i,j,k)-t(i,j-1,k);
        end
    end
end
% define the degradation rates
v = v;
% define the covariance matrix
sigma2 = sigma2; 
% define parameters in error distribution
eta = eta;
for i = 1:n
    for j = 1:m
        zeta(i,j) = gamrnd(eta,1./eta); 
        % zeta(i,j) is varsigma_{ij}
    end
end
% let re-think the model and simulate more accurate
for i = 1:1:n
    for j = 2:1:m
        for k = 1:1:K
            dY(i,1,k,1:d) = zeta(i,1)*mvnrnd(v*dt(i,1,k),sigma2*dt(i,1,k))';
            Y(i,1,k,1:d) = dY(i,1,k,1:d);
            % first simulate the degradation increment, then calculate the
            % degradatoin observation
            dY(i,j,k,1:d) = zeta(i,j)*mvnrnd(v*dt(i,j,k),sigma2*dt(i,j,k))';
            Y(i,j,k,1:d) = Y(i,j-1,k,1:d)*zeta(i,j)/zeta(i,j-1)+dY(i,j,k,1:d);
        end
    end
end
if debug_flag == 1
    plot(t,Y);
else
end
end