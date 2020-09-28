% Coarsened posterior distribtuion of covariance matrix (Sigma) with Jefferys prior as
% non-informative prior of Sigma
% Note: (1) The discrepancy between the empirical desnity of oberserved and the empirical
%           density of idealized samples is quantified using Kullback-Liebler
%           divergence
%       (2) The power posterior approximation (valid only for small
%           discrepancies) is used
%       (3) the coarsening parameter \alpha (Miller and Dunson, 2019)
%           should be such that degrees of freedom of of inverse-Wishart
%           distribution is greater than dimensionality of the data
% Ref: Miller and Dunson (2019). Robust Bayesian inference via coarsening.

clear all
close all
clc

direc='D:/Research/Thesis_work/Structural_vs_measurement_uncertainty/matlab_codes';
N=1000; % number of covariance matrices to be drawn
alpha=1; % coarsening parameter


% load observed streamflow data (ML estimated)
fname='streamflow_18.mat';
filename=fullfile(direc,'huc_04100003','SJRW_ML_results',fname);
load(filename);                     % streamflow in cfs
strm_ML=MC1_temp*0.028316847;       % streamflow in cms
strm_datenums = datenum('01-01-2005','dd-mm-yyyy'):datenum('31-12-2016','dd-mm-yyyy');
begin_datenum = datenum('01-01-2008','dd-mm-yyyy');
end_datenum = datenum('31-12-2010','dd-mm-yyyy');
ind = find(strm_datenums>=begin_datenum & strm_datenums<=end_datenum);
strm_ML=strm_ML(ind,:);            % retain only one year of data

% define parameters of prior  (inverse Wishart distribtuion)
[d,n]=size(strm_ML); % d=dimensions, n=number of observations
nu=d+1;           % degrees of freedom
Lambda=1*eye(d);        % scale matrix

% compute \zeta_{n} parameter (Miller and Dunson, 2019)
zeta_n=alpha/(alpha+n);

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

% define parameters of posterior inverse wishart distribution
nu_s=nu+zeta_n*n;
Lambda_s=Lambda+zeta_n*Q;

% draw samples from inverse Wishart distribtuion
C = zeros(d,d,N);
for samp_ind=1:N
    
    C(:,:,samp_ind)=iwishrnd(Lambda_s,nu_s);         % samples of covariance matrix from posterior distribution
%     Cr(:,:,samp_ind)=corrcov(C(:,:,samp_ind));       % samples of correlation matrix from posterior distribution
       
end

% save the samp and R
fname=strcat('covmat','_','alpha = ',num2str(alpha),'.mat');
save_filename=fullfile('E:/covmat_Jefferys_prior',fname);
save(save_filename,'C','-v7.3');