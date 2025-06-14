function [v1_ini,v2_ini,v3_ini,v_ini,sigma2_ini,eta_ini,zeta_ini] = Wiener_MCEM_initiation_3d(dY,dt,n,K,m)
%{ 

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

- Model:
- use basic multivariate Wiener process model to initialize v and sigma2
- specifically: X_{ik}(t) = vt+ sigma2^(1/2)* B(t)


- Notations:
- i is the index for the i-th batch
- k is the index for the k-th unit in batch i
- j is the index for the j-th measurement

%} 

for i = 1:n
    for k = 1:K
        for j = 2:m
            v_ini_d1(i,j,k) = dY(i,j,k,1)./dt(i,j,k);
            v_ini_d2(i,j,k) = dY(i,j,k,2)./dt(i,j,k);
            v_ini_d3(i,j,k) = dY(i,j,k,3)./dt(i,j,k);
        end
    end
end
v1_ini = sum(sum(sum(v_ini_d1)))/(n*K*(m-1));
v2_ini = sum(sum(sum(v_ini_d2)))/(n*K*(m-1));
v3_ini = sum(sum(sum(v_ini_d3)))/(n*K*(m-1));
v_ini = [v1_ini,v2_ini,v3_ini]';

sigma2_ini = [0 0 0; 0 0 0; 0 0 0];

for i = 1:n
    for j = 2:m
        for k = 1:K
            dy = [dY(i,j,k,1);dY(i,j,k,2);dY(i,j,k,3)];
            sigma2_ini = sigma2_ini + (dy-v_ini*dt(i,j,k))*(dy-v_ini*dt(i,j,k))'/dt(i,j,k);
        end
    end
end

sigma2_ini = sigma2_ini/(n*(m-1)*K);


% find the initial value of eta
for i = 1:n
    for j = 2:m
        zeta_ini(i,j) = (mean(v_ini_d1(i,j,:))/v1_ini+...
            mean(v_ini_d2(i,j,:))/v2_ini+...
            mean(v_ini_d3(i,j,:))/v3_ini)/3;
    end
end
varsigma_now = reshape(zeta_ini,[],1);
varsigma_ini_nozero = varsigma_now(varsigma_now>0);
eta_ini = length(varsigma_ini_nozero)*sum(varsigma_ini_nozero)./...
    (length(varsigma_ini_nozero)*sum(varsigma_ini_nozero.*log(varsigma_ini_nozero))-sum(varsigma_ini_nozero).*sum(log(varsigma_ini_nozero)));
end

