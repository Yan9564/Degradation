function [opts] = Wiener_MCEM_3d(dt,dY,n,m,d,MCEM_set_params,K,v_ini,sigma2_ini,eta_ini,zeta_ini)
%{

- Description:
- this function estimates the model parameters in Wiener process with 3 PCs
using MCEM algorithm.
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
- opts: a struct saving all outputs
- - opts.eta_est: estimates of eta
- - opts.v_est: estimates of v, i.e., the degradation rates
- - opts.v1_est: estimates of v1, i.e., degradation rate of the first PC
- - opts.v2_est: estimates of v2, i.e., degradation rate of the second PC
- - opts.v3_est: estimates of v3, i.e., degradation rate of the third PC
- - opts.sigma2_est: estimates of sigma
- - opts.sigma2_11_est: estimates of sigma_11
- - opts.sigma2_12_est: estimates of sigma_12
- - opts.sigma2_13_est: estimates of sigma_13
- - opts.sigma2_22_est: estimates of sigma_22
- - opts.sigma2_23_est: estimates of sigma_23
- - opts.sigma2_33_est: estimates of sigma_33
- - opts.zeta_est: estimates of varsigma, i.e., the measurement error

%}

% strip off the parameters in struct
Names = fieldnames(MCEM_set_params);% List all variables under MyStructure
for nn = 1:length(Names)
    eval([Names{nn},' = MCEM_set_params.',Names{nn},';']);% Assign data to original names
end

% the MCEM algorithm for estimating zeta v Σ
% initialization for v
v1_ini = v_ini(1);
v2_ini = v_ini(2);
v3_ini = v_ini(3);
v_est{1} = v_ini;
% initialization for sigma
sigma2_ini;
sigma2_est{1} = sigma2_ini;
% initialization for eta
eta_ini;
eta_est_ex(1) = eta_ini;
% initialization for zeta
zeta_ini;
for l = 1
    for i = 1:n
        for j = 1:m
            for b = 1
                zeta_est(l,i,j,b) = zeta_ini(i,j);
                zeta_now(l,i,j,b) = zeta_ini(i,j);
            end
        end
    end
end

% number of Gibbs sampling
B;
% number of MCEM iteration
L;
% tolerance of MCEM iteration
tol;

for l = 1:L
    % E-step
    % MCMC Gibbs sampling
    fprintf('===== begin the Wiener_MCEM at iteration:')
    disp(l)
    eta_est(1) = eta_ini;
    for i = 1:n
        for j = 2:m
            count = 0;
            zeta_acpt = [];
            zeta_rjct = [];
            for b = 2:10000
                % initialization
                %{
                    
                    zeta_now(l,i,1,2):
                    the varisigma at the l-th iteration of MCEM
                    for the i-th batch
                    at the 1-st measurement
                    in the 2-nd Gibbs sampling
                    
                    zeta_est is the Gibbs samples of zeta;
                    zeta_now is the Gibbs samples of zeta after thinned operation;
                    zeta_new is the potential candidate of zeta
                    
                %}
                
                mcmc_sigma_c = mcmc_sigma_c;
                zeta_ran = normrnd(0,mcmc_sigma_c.^2);
                zeta_new(l,i,j,b) = zeta_est(l,i,j,b-1)+zeta_ran;
                
                if zeta_new(l,i,j,b)<=0
                    zeta_new(l,i,j,b) = zeta_est(l,i,j,b-1);
                    % to make sure that zeta is in the positive support
                else
                end
                % lnf_new = lnf1_new+lnf2_new
                % lnf1_new is the likelihood of zeta
                % lnf2_new is the likelihood of dy
                
                lnf1_new(l,i,j,b) = (eta_est_ex(l)-1)*log(zeta_new(l,i,j,b))...
                    - eta_est_ex(l)*zeta_new(l,i,j,b);
                lnf1_old(l,i,j,b) = (eta_est_ex(l)-1)*log(zeta_est(l,i,j,b-1))...
                    - eta_est_ex(l)*zeta_est(l,i,j,b-1);
                
                for k = 1:K
                    dy = [dY(i,j,k,1);dY(i,j,k,2);dY(i,j,k,3)];
                    
                    lnf2_new(l,i,j,k,b) = -0.5*...
                        (dy/zeta_new(l,i,j,b)-v_est{l}*dt(i,j,k))'*...
                        pinv(sigma2_est{l}*dt(i,j,k))*...
                        (dy/zeta_new(l,i,j,b)-v_est{l}*dt(i,j,k));
                    
                    lnf2_old(l,i,j,k,b) = -0.5*...
                        (dy/zeta_est(l,i,j,b-1)-v_est{l}*dt(i,j,k))'*...
                        pinv(sigma2_est{l}*dt(i,j,k))*...
                        (dy/zeta_est(l,i,j,b-1)-v_est{l}*dt(i,j,k));
                end
                
                % calculate the acceptance probability
                
                f_new(l,i,j,b) = exp(lnf1_new(l,i,j,b)+sum(lnf2_new(l,i,j,:,b)));
                f_old(l,i,j,b) = exp(lnf1_old(l,i,j,b)+sum(lnf2_old(l,i,j,:,b)));
                h(l,i,j,b) = f_new(l,i,j,b)/f_old(l,i,j,b);
                
                r = 0.2;
                
                if r <= h(l,i,j,b)
                    zeta_acpt = [zeta_acpt, zeta_new(l,i,j,b)];
                    zeta_est(l,i,j,b) = zeta_new(l,i,j,b);
                    count = count+1;
                    % fprintf('accepted new zeta \n')
                else
                    zeta_rjct = [zeta_rjct, zeta_est(l,i,j,b-1)];
                    zeta_est(l,i,j,b) = zeta_est(l,i,j,b-1);
                    % fprintf('rejected new zeta \n')
                end

                if count/1000>0.5
                    break
                end
            end
            accept_rate(l,i,j) = count/B;
            % thinned MCMC sequence
            rate = 0.45;
            start = floor(mcmc_burn_percent*B);
            start_acpt = floor(start*rate);
            start_rjct = floor(start*(1-rate));
            M = length([start:mcmc_interval:B]);
            M_acpt = floor(rate*M);
            M_rjct = M-M_acpt;
            mcmc_inte_acpt = floor(mcmc_interval*rate);
            mcmc_inte_rjct = mcmc_interval - mcmc_inte_acpt;
            zeta_now_sequence = [zeta_rjct(start_rjct:mcmc_inte_rjct:start_rjct+M_rjct*mcmc_inte_rjct),zeta_acpt(start_acpt:mcmc_inte_acpt:start_acpt+M_acpt*mcmc_inte_acpt)];
            zeta_now(l,i,j,1:M) = zeta_now_sequence(1:M);
            zeta_est_EX(l,i,j) = mean(zeta_now(l,i,j,:));
            % initialization of zeta for (l+1)-th iteration of MCEM
            zeta_est(l+1,i,j,1) = zeta_now(l,i,j,M);
        end
    end
    
    % M-step
    % the actual degradation proecss X
    for i = 1:n
        for j = 2:m
            for k = 1:K
                dX(l,i,j,k,1,1:M) = dY(i,j,k,1)./zeta_now(l,i,j,1:M);
                dX(l,i,j,k,2,1:M) = dY(i,j,k,2)./zeta_now(l,i,j,1:M);
                dX(l,i,j,k,3,1:M) = dY(i,j,k,3)./zeta_now(l,i,j,1:M);
            end
        end
    end
    % calculate the drift parameter
    v1_est(l+1) = sum(sum(sum(sum(dX(l,:,2:m,:,1,:)))))/...
        (M*sum(sum(sum(dt(:,2:m,:)))));
    v2_est(l+1) = sum(sum(sum(sum(dX(l,:,2:m,:,2,:)))))/...
        (M*sum(sum(sum(dt(:,2:m,:)))));
    v3_est(l+1) = sum(sum(sum(sum(dX(l,:,2:m,:,3,:)))))/...
        (M*sum(sum(sum(dt(:,2:m,:)))));
    v_est{l+1} = [v1_est(l+1);v2_est(l+1);v3_est(l+1)];
    % estimate the sigma2
    sigma2_est{l+1} = [0,0,0;0,0,0;0,0,0];
    for i = 1:n
        for j = 2:m
            for k = 1:K
                for b = 1:M
                    dX_vec = [dX(l,i,j,k,1,b);
                        dX(l,i,j,k,2,b);
                        dX(l,i,j,k,3,b)];
                    sigma2_est{l+1}(1:d,1:d) = sigma2_est{l+1}(1:d,1:d) +...
                        (dt(i,j,k))^(-1).*...
                        (dX_vec-v_est{l+1}*dt(i,j,k))*...
                        (dX_vec-v_est{l+1}*dt(i,j,k))';
                end
            end
            zeta_ex(l,i,j) = zeta_now(l,i,j,M);
        end
    end
    sigma2_est{l+1}(1:d,1:d) = sigma2_est{l+1}(1:d,1:d)/(M*K*(m-1)*n);
    sigma2_11_est(l+1) = sigma2_est{l+1}(1,1);
    sigma2_12_est(l+1) = sigma2_est{l+1}(1,2);
    sigma2_13_est(l+1) = sigma2_est{l+1}(1,3);
    sigma2_22_est(l+1) = sigma2_est{l+1}(2,2);
    sigma2_23_est(l+1) = sigma2_est{l+1}(2,3);
    sigma2_33_est(l+1) = sigma2_est{l+1}(3,3);
    
    % estimate eta
    
    zeta_eta = reshape(zeta_now(l,:,:,:),[],1);
    zeta_eta = zeta_eta(zeta_eta>0);
    M_zeta = length(zeta_eta);
    eta_est(l+1) = M_zeta*sum(zeta_eta)/...
        (M_zeta*sum(zeta_eta.*log(zeta_eta))-...
        sum(log(zeta_eta))*sum(zeta_eta));
    
    % calaulte eta expectation
    zeta_ex_eta = reshape(zeta_ex(l,:,:),[],1);
    zeta_ex_eta = zeta_ex_eta(zeta_ex_eta>0);
    M_zeta_ex = length(zeta_ex_eta);
    eta_est_ex(l+1) = M_zeta_ex*sum(zeta_ex_eta)/...
        (M_zeta_ex*sum(zeta_ex_eta.*log(zeta_ex_eta))-...
        sum(log(zeta_ex_eta))*sum(zeta_ex_eta));
    
    
    % check the difference between two iterations
    v1_est_diff(l+1) = v1_est(l+1)-v1_est(l);
    v2_est_diff(l+1) = v2_est(l+1)-v2_est(l);
    v3_est_diff(l+1) = v3_est(l+1)-v3_est(l);
    sigma2_11_est_diff(l+1) = sigma2_est{l+1}(1,1)-sigma2_est{l}(1,1);
    sigma2_22_est_diff(l+1) = sigma2_est{l+1}(2,2)-sigma2_est{l}(2,2);
    sigma2_33_est_diff(l+1) = sigma2_est{l+1}(3,3)-sigma2_est{l}(3,3);
    sigma2_12_est_diff(l+1) = sigma2_est{l+1}(1,2)-sigma2_est{l}(1,2);
    sigma2_13_est_diff(l+1) = sigma2_est{l+1}(1,3)-sigma2_est{l}(1,3);
    sigma2_23_est_diff(l+1) = sigma2_est{l+1}(2,3)-sigma2_est{l}(2,3);
    eta_est_diff(l+1) = eta_est_ex(l+1)-eta_est_ex(l);
    % check the tolerance
    if abs(v1_est_diff(l+1))<tol...
            && abs(v2_est_diff(l+1))<tol ...
            && abs(v3_est_diff(l+1))<tol ...
            && abs(sigma2_11_est_diff(l+1))<tol...
            && abs(sigma2_22_est_diff(l+1))<tol...
            && abs(sigma2_33_est_diff(l+1))<tol...
            && abs(sigma2_12_est_diff(l+1))<tol...
            && abs(sigma2_13_est_diff(l+1))<tol...
            && abs(sigma2_23_est_diff(l+1))<tol...
            && abs(eta_est_diff(l+1))<tol
        mcem_flag = l;
        fprintf('tolerance satisfied \n')
        break;
    else
        mcem_flag = 0;
    end
end
opts = struct();
opts.v_est = v_est;
opts.sigma2_est = sigma2_est;
opts.eta_est = eta_est;
opts.zeta_est = zeta_est;
opts.zeta_now = zeta_now;
opts.v1_est = v1_est;
opts.v2_est = v2_est;
opts.v3_est = v3_est;
opts.sigma2_11_est = sigma2_11_est;
opts.sigma2_12_est = sigma2_12_est;
opts.sigma2_13_est = sigma2_13_est;
opts.sigma2_22_est = sigma2_22_est;
opts.sigma2_23_est = sigma2_23_est;
opts.sigma2_33_est = sigma2_33_est;
opts.mcem_flag = mcem_flag;
end


