% this routine draws samples from distribution obtained by random forest streamflow
% time-series 

clear all
close all
clc

direc = 'D:/Research/Thesis_work/Structural_vs_measurement_uncertainty/matlab_codes';
alpha = 1;


% read streamflow data
fname = 'streamflow_18.mat';
filename = fullfile(direc,'huc_04100003/SJRW_ML_results',fname);
load(filename);
MC1_temp = MC1_temp*0.0283; % cfs to cms
strm_ML = mean(MC1_temp,2);
strm_datenums = datenum('01-01-2005','dd-mm-yyyy'):datenum('31-12-2016','dd-mm-yyyy');
begin_datenum = datenum('01-01-2008','dd-mm-yyyy');
end_datenum = datenum('31-12-2010','dd-mm-yyyy');
ind = find(strm_datenums>=begin_datenum & strm_datenums<=end_datenum);
strm_ML = strm_ML(ind);
strm_datenums = strm_datenums(ind);

% read covraince matrix data
fname = ['covmat_alpha =',num2str(alpha),'.mat'];
filename = fullfile('E:/covmat_Jefferys_prior',fname);
load(filename);

% draw streamflow samples for each covaraince matrix
dir_name = ['strm_samples_alpha=',num2str(alpha)];
mkdir(fullfile('E:/covmat_Jefferys_prior',dir_name))
for ci = 1:size(C,3)
    
    Ctmp = C(:,:,ci);
    strm_tmp = mvnrnd(strm_ML,Ctmp,1000);
    strm_tmp(strm_tmp<0) = 0;
    sname = ['strm_samples_alpha=',num2str(alpha),'_',num2str(ci),'.mat'];
    filename = fullfile('E:/covmat_Jefferys_prior',dir_name,sname);
    save(filename,'strm_tmp','-v7.3');
end

