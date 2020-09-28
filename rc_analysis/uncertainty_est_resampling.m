% this routine resamples 100000 samples from streamflow population obtained
% by rating curve analysis

clear all
close all
clc

direc = 'D:/Research/Thesis_work/Structural_vs_measurement_uncertainty/matlab_codes';

% read data
fname = 'rc_strm_samples_04180500.mat';
filename = fullfile(direc,'huc_04100003/results',fname);
load(filename)

% draw samples
probs = freq'/sum(freq);
samp_space = (1:length(probs))';
strm_samps = zeros(100000,size(strm_samples,1));

for count = 1:100000
    
    % sample covariance matrix
    r = prob_discrete(samp_space,probs);
    % sample a row corresponding to a covariance matrix
    rind = randsample(1:1000,1);
    strm_samps(count,:) = strm_samples(:,rind,r)';

end

% save data
sname = 'rc_strm_samples_04180500_resampled.mat';
filename = fullfile(direc,'huc_04100003/results',sname);
save(filename,'strm_samps');
