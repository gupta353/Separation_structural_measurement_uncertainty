% Posterior distribtuion of covariance matrix (Sigma) with Jefferys prior as
% non-informative prior of Sigma

% clear all
% close all
% clc

direc='D:/Research/Thesis_work/Structural_vs_measurement_uncertainty/matlab_codes';
N=1000; % number of covariance matrices to be drawn

% load observed streamflow data (ML estimated)
fname='streamflow_200.mat';
filename=fullfile(direc,'huc_04100003','SJRW_ML_results',fname);
load(filename);                     % streamflow in cfs
strm_ML=MC1_temp*0.028316847;       % streamflow in cms
strm_ML=strm_ML(1:365,:);            % retain only one year of data
break
% define parameters of prior  (inverse Wishart distribtuion)
[d,n]=size(strm_ML); % d=dimensions, n=number of observations
nu=0;           % degrees of freedom
Lambda=zeros*eye(d);        % scale matrix

% compute error in each realization
mean_strm_ML=mean(strm_ML,2);
error_ML=bsxfun(@minus,strm_ML,mean_strm_ML);


%% compute posterior distribtuion
% compute of scale matrix
Q=zeros(d,d);       % scale matrix
for Q_ind=1:n
    
    Q=Q+error_ML(:,Q_ind)*error_ML(:,Q_ind)';
    
end

% compute determinant of Q
eigenvalues=eig(Q);
log_deter_Q=sum(log(abs(eigenvalues)));
deter_Q=exp(log_deter_Q);
break
% define parameters of posterior inverse wishart distribution
nu_s=nu+n;
Lambda_s=Lambda+Q;

% draw samples from inverse Wishart distribtuion
for samp_ind=1:N
    
    C(:,:,samp_ind)=iwishrnd(Lambda_s,nu_s);         % samples of covariance matrix from posterior distribution
    Cr(:,:,samp_ind)=corrcov(C(:,:,samp_ind));        % samples of correlation matrix from posterior distribution
       
end

% save the samp and R
fname='covmat.mat';
save_filename=fullfile(direc,'huc_04100003','results','covmat_Jefferys_prior',fname);
save(save_filename,'C','Cr');