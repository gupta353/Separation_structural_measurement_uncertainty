% computation of uncertainty bounds over streamflow measurement errors
% using rating curve analysis

clear all
close all
clc

direc='D:/Research/Thesis_work/Structural_vs_measurement_uncertainty/matlab_codes';
begin_date='2001-01-01';    % begin date in yyyy-mm-dd format
end_date='2013-12-31';      % end date in yyyy-mm-dd format

% read the .mat file containing MCMC chains
fname = 'rc_MCMC_chain_04180500.mat';
filename = fullfile(direc,'huc_04100003/results',fname);
load(filename);

beg_samp = 12001;               % begin index of sample such that Rd is less than 1.1
nchains = size(smpl,3);         % number of chains
smpl = smpl(beg_samp:end,:,:);
smpl_tmp=smpl;
smpl=[];
% mix all the samples from different chains
for ind = 1:nchains
    smpl = [smpl;smpl_tmp(:,:,ind)];
end

% select 10000 samples randomly out of smpl
nsamps = size(smpl,1);
inds = randsample(1:nsamps,10000);
smpl = smpl(inds,:);

% read gage height data
fname = 'gage_height_04180500.txt';
filename = fullfile(direc,'huc_04100003/streamflow_data',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);
gh_dates = data{1};
gage_height = data{2}*0.3048;
gh_datenums = cellfun(@(x)datenum(x,'mm/dd/yyyy HH:MM'),gh_dates);

% read streamflow data
fname = 'streamflow_04180500.txt';
filename = fullfile(direc,'huc_04100003/streamflow_data',fname);
fid = fopen(filename,'r');
data = textscan(fid,'%s%f','delimiter','\t','headerlines',1);
fclose(fid);
strm_dates = data{1};
strm = data{2}*0.0283;

begin_datenum = datenum('01-01-2008','dd-mm-yyyy');
end_datenum = datenum('31-12-2010','dd-mm-yyyy');
strm_datenums = cellfun(@(x)datenum(x,'yyyy-mm-dd'),strm_dates);
ind = find(strm_datenums>=begin_datenum & strm_datenums<=end_datenum);
strm = strm(ind);
strm_datenums = strm_datenums(ind);

% compute streamflow uncertianty for each parameter set in smpl
uni_smpl=unique(smpl,'rows','stable');
aspace = 'log';
strm_samples = zeros(length(strm_datenums),1000,size(uni_smpl,1));
for theta_ind=1:size(uni_smpl,1)
    theta=uni_smpl(theta_ind,:);
    m       =    theta(1);                     % number of segments
    h01     =    theta(2);                     % cease-to-flow parameter of the first segment
    h_s     =    theta(2:m+2);                 % list of break points (h01 is included in the list of break-points)
    h0_list =    [h01,theta(m+3:2*m+1)];       % list of cease-to-flow parameters
    a1      =    theta(2*m+2);                 % list of multiplier parameters
    b_list  =    theta(2*m+3:3*m+2);           % list of exponent parameters
    sigma2_list  =    theta(3*m+3:4*m+2);
    
    log_Q_sim = rc_est(gage_height,h_s,a1,b_list,h0_list,aspace);
    Q_sim = exp(log_Q_sim);
    
    % assign standard deviations to Q_sim
    for qsim_ind = 1:length(Q_sim)
        seg = find(gage_height(qsim_ind)<h_s,1,'first')-1;
        var_Q_sim(qsim_ind,1) = sigma2_list(seg);
    end
    
    
    % compute average daily streamflow
    for dind = 1:length(strm_datenums)
        
        ind = find(gh_datenums>=strm_datenums(dind) & gh_datenums<=strm_datenums(dind)+1);
        Qd_sim(dind) = mean(exp(log_Q_sim(ind)));
        var_Qd_sim(dind) = sum(var_Q_sim(ind))/length(ind)^2;
        strm_samples(dind,:,theta_ind) = normrnd_truncated(Qd_sim(dind),sqrt(var_Qd_sim(dind)),1000)';
        
    end
    
    % determine the frequency of current parameter set
    freq(theta_ind) = sum(ismember(smpl,theta,'rows'));
    
end
save('rc_strm_samples_04180500.mat','freq','strm_samples','-v7.3');