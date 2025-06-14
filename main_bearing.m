%{

- Description: 
- This file analyzes the degradation data from bearing 1 and bearing 2.
- See details of the analysis in supplementary in "Bearings on high-speed trains".

- Input:
- bearing_1.mat and bearing_2.mat: degradation data on Bearing 1 and Bearing 2 from high-speed train.

- Output:
- parameter_estimation_proposed_model_bearing.mat: parameter estimates at each iteration of the MCEM algorithm.

- Notation:
- in the code  | in the paper    | Explanations
- d            | d               | number of PCs in each unit.
- B            | B_0+Mu          | number of MCMC samples in Algorithm 1.
- n            | n               | number of batch size in the degradation test.
- MCEM_N       | \               | number of iterations in MCEM algorithm.
- m            | m               | the total number of measurements in batch i.
- Y            | Y               | degradation levels of the d PCs. 
- t            | t               | measurement times.
- M0           | \               | basic multivariate Wiener process model without multiplicative measurement errors.
- v            | v               | degradation rates in Wiener process.
- sigma2       | Sigma           | covariance matrix in Wiener process.
- eta          | eta             | both the shape and rate parameters in gamma distribution.
- _est         | \               | estimation of the parameters.


%}


clc;
clear all;
close all;

% load degradation data
d = 2;
MCEM_N = 100;
B = 1000;
n = 2; 
Bearing = 1; % index of the bearing
filename = ['bearing_' num2str(Bearing) '.mat'];
load(filename);
db(1,:) = db_value_feature(3:22);
sv(1,:) = sv_value_feature(3:22);
Bearing = 2;
filename = ['bearing_' num2str(Bearing) '.mat'];
load(filename);
db(2,:) = db_value_feature(3:22);
sv(2,:) = sv_value_feature(3:22);
% take the logarithm of the SPM-based shock pulse signals
log_db = log(db);
log_sv = log(sv);
t=time_feature(3:22);
% smooth the degradation path
log_db(1,:) = smooth(log_db(1,:),5);
log_sv(1,:) = smooth(log_sv(1,:),5);
log_db(2,:) = smooth(log_db(2,:),5);
log_sv(2,:) = smooth(log_sv(2,:),5);
% prepare the data for estimation
m = length(t);
for i=1:n
    for j=1:m
        Y{i,j}=[log_db(i,j),log_sv(i,j)]';
    end
end
for i=1:n
    DY{i,1}=[0,0]';
    for j=2:m
        DY{i,j}=Y{i,j}-Y{i,j-1};
        Dt(i,j)=t(j)-t(j-1);
    end
end
for i=1:n
    for j=1:m
        Y_now1(j)=Y{i,j}(1,1);
        Y_now2(j)=Y{i,j}(2,1);
    end
end
% estimate model parameters using basic multivariate Wiener process model 
for i=1:n
    for j=2:m
        DY_M0(:,j+m*(i-1))=DY{i,j};
        DY1_M0(:,j+m*(i-1))=DY{i,j}(1,1);
        DY2_M0(:,j+m*(i-1))=DY{i,j}(2,1);
        dt_M0(j+m*(i-1))=t(j)-t(j-1);
    end
end
v_M0(:,1)=sum(DY_M0,2)./sum(dt_M0);
for i=1:n
    for j=2:m
        sigma2_M0_sub(:,:,j+(i-1)*m)=(DY_M0(:,j+(i-1)*m)-v_M0*(t(j)-t(j-1)))*...
            (DY_M0(:,j+(i-1)*m)-v_M0*(t(j)-t(j-1)))'/(n*(m-1)*(t(j)-t(j-1)));
    end
end
for k=1:d
    for p=1:d
        sigma2_M0_est(k,p)=sum(sigma2_M0_sub(k,p,:));
    end
end
% Estimate the model parameters using MCEM
eta_est(1) = 100;
v_est{1} = v_M0;
sigma2_est{1} = sigma2_M0_est;
for l=1:MCEM_N
    % E-step
    % MCMC sampling
    for i=1:n
        for j=2:m
            count=0;
            for b=2:B
                % set initial values for zeta
                zeta_now(l,i,1,2)=1;  
                zeta_est_EX(l,i,1)=1;  
                zeta_est(1,1,2,1)=1;  
                zeta_est_EX_inv(l,i,1)=1;  
                
                sigma_c = sqrt(0.1); % e tuning parameter sigma_c in Algorithm 1 in the paper
                zeta_ran=normrnd(0,sigma_c^2);
                zeta_new(l,i,j,b)=zeta_est(l,i,j,b-1)+zeta_ran;  
                
                if zeta_new(l,i,j,b)<=0
                    zeta_new(l,i,j,b)=zeta_est(l,i,j,b-1); % make sure that zeta is larger than 0
                else
                end
                r = 0.98; % tuning parameter r in MCMC Algorithm 1 
                % calculate value of h function in Algorithm 1 in the paper
                h_new(l,i,j,b)=(eta_est(l)-1)*log(zeta_new(l,i,j,b))-eta_est(l)*zeta_new(l,i,j,b)-0.5*...
                    (DY{i,j}/zeta_new(l,i,j,b)-v_est{l}*(t(j)-t(j-1)))'*...
                    inv(sigma2_est{l}*(t(j)-t(j-1)))*...
                    (DY{i,j}/zeta_new(l,i,j,b)-v_est{l}*(t(j)-t(j-1)));
                
                h_old(l,i,j,b)=(eta_est(l)-1)*log(zeta_est(l,i,j,b-1))-eta_est(l)*zeta_est(l,i,j,b-1)-0.5*...
                    (DY{i,j}/zeta_est(l,i,j,b-1)-v_est{l}*(t(j)-t(j-1)))'*...
                    inv(sigma2_est{l}*(t(j)-t(j-1)))*...
                    (DY{i,j}/zeta_est(l,i,j,b-1)-v_est{l}*(t(j)-t(j-1)));
                
                beta(l,i,j,b)=h_old(l,i,j,b)/h_new(l,i,j,b);
                if isnan(beta(l,i,j,b)) == 1
                    alpha(l,i,j,b)=0;
                else
                    alpha(l,i,j,b)=min(beta(l,i,j,b),1);
                end
                
                if r<=alpha(l,i,j,b)
                    zeta_est(l,i,j,b)=zeta_new(l,i,j,b);
                    count=count+1;
                else
                    zeta_est(l,i,j,b)=zeta_est(l,i,j,b-1);
                end
            end
            accept_rate(l,i,j)=count/B;
            
            % thinned MCMC sequence
            start=floor(0.2*B);
            interval=10;
            q=1;
            for g=start:interval:B
                zeta_now(l,i,j,q)=zeta_est(l,i,j,g);
                q=q+1;
            end
            
            zeta_est_EX(l,i,j)=mean(zeta_now(l,i,j,:));
            zeta_est_EX_inv(l,i,j)=mean(zeta_now(l,i,j,:).^(-1));
            zeta_est_EX2(l,i,j)=mean(zeta_now(l,i,j,:).^2);
            zeta_now_size(l,i,j)=length(zeta_now(l,i,j,:));
            zeta_est(l,i,j+1,1)=zeta_est_EX(l,i,j);
        end
        zeta_est(l,i+1,2,1)=zeta_est_EX(l,i,2);
    end
    zeta_est(l+1,1,2,1)=zeta_est_EX(l,1,2);
    
    % M-step
    for i=1:n
        for j=2:m
            DY_now(:,j+(i-1)*m)=DY{i,j};
            zeta_est_EX_inv_now(j+(i-1)*m)=zeta_est_EX_inv(l,i,j);
            zeta_est_EX_inv_now(1+(i-1)*m)=1;
            zeta_est_EX_now(j+(i-1)*m)=zeta_est_EX(l,i,j);
            zeta_est_EX_now(1+(i-1)*m)=1;
            zeta_est_EX2_now(j+(i-1)*m)=zeta_est_EX2(l,i,j);
            dt_now(j+(i-1)*m)=t(j)-t(j-1);
            DX_now(:,j+(i-1)*m)=DY{i,j}*zeta_est_EX_now(j+(i-1)*m);
        end
    end
    v_est{l+1}=sum(DX_now,2)/sum(dt_now);
    
    for i=1:n
        for j=2:m
            for b=1:zeta_now_size(l,i,j)
                zeta_now(l,i,1,b)=1;
                sigma2_sub(l,:,:,j+(i-1)*m,b)=(DY{i,j}/zeta_now(l,i,j,b)-v_est{l+1}*dt_now(j+(i-1)*m))*...
                    (DY{i,j}/zeta_now(l,i,j,b)-v_est{l+1}*dt_now(j+(i-1)*m))'/...
                    dt_now(j+(i-1)*m);
            end
        end
    end
    for k=1:d
        for p=1:d
            sigma2_est{l+1}(k,p)=sum(sum(sigma2_sub(l,k,p,2:end,:)))./sum(sum(zeta_now_size(l,:,2:end)));
        end
    end
    
    zetas=reshape(zeta_now(l,:,:,:),[],1);
    zetas_nozero=zetas(zetas>0);
    % estimate eta in the gamma distirbution 
    eta_est(l+1)=length(zetas_nozero)*sum(zetas_nozero)./(length(zetas_nozero)*sum(zetas_nozero.*log(zetas_nozero))-sum(zetas_nozero).*sum(log(zetas_nozero)));
end

% Parameter estimates at each iteration of the MCEM algorithm
for l=1:MCEM_N
    v1_est(l)=v_est{l}(1);
    v2_est(l)=v_est{l}(2);
    sigma2_11_est(l)=sigma2_est{l}(1,1);
    sigma2_12_est(l)=sigma2_est{l}(1,2);
    sigma2_22_est(l)=sigma2_est{l}(2,2);
end
v_result(:) = v_est{MCEM_N};
sigma2_result(:,:) = sigma2_est{MCEM_N};
eta_result(:) = eta_est(MCEM_N);

% calculate AIC values for the proposed model fitting the bearings.
for i=1:n
    for j=2:m
        zeta(i,j)=mean(zeta_now(MCEM_N,i,j,:));
        L_sub(i,j)= 0.5*log(det(sigma2_result*Dt(i,j)))+...
            0.5*(DY_M0(:,j+(i-1)*m)-v_result'*Dt(i,j))'*inv(sigma2_result*Dt(i,j))*(DY_M0(:,j+(i-1)*m)-v_result'*Dt(i,j))+...
            (eta_result.*log(eta_result)-log(gamma(eta_result))+(eta_result-1).*log(zeta(i,j))-eta_result.*zeta(i,j));
    end
end
% the result of L is shown in "Table: AIC values for candidate models
% fitting the bearings" in the Section of "Supplementary Results"
L = sum(sum(L_sub));

% save the Parameter estimates at each iteration of the MCEM algorithm
filename=['results\parameter_estimation_proposed_model_bearing.mat'];
save(filename,'t','v_result','sigma2_result','eta_result','v1_est','v2_est','sigma2_11_est','sigma2_12_est','sigma2_22_est','eta_est','L');
