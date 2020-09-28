% Visualization of covariance matrix distribution

clear all
close all
clc

% time-steps at which boxplot is to plotted
boxplot_t=1:10:365;

% load data
direc='D:/Research/Thesis_work/Structural_vs_measurement_uncertainty/matlab_codes';
fname='covmat.mat';
filename=fullfile(direc,'huc_04100003','results','covmat_Jefferys_prior',fname);
load(filename);         % covariance matrix and correlation matrices


% plot boxplots of standard deviations
for sig_ind=1:size(C,1)
    
    sigma(:,sig_ind)=sqrt(squeeze(C(sig_ind,sig_ind,:)));
        
end
figure;
boxplot(sigma(:,boxplot_t),boxplot_t);
xlabel('Time-step','fontname','arial','fontsize',12);
ylabel('Standard deviation (m^3 s^{-1})','fontname','arial','fontsize',12);
set(gca,'fontname','arial','fontsize',12);

% plot boxplots of correlations
for corr_ind=2:size(C,1)
    
    corrlat(:,corr_ind)=squeeze(Cr(1,corr_ind,:));
        
end
figure;
boxplot(corrlat(:,boxplot_t+1),boxplot_t+1);
xlabel('n^{th} time-step','fontname','arial','fontsize',12);
ylabel('Correlation coefficient between 1^{st} and n^{th} time-step','fontname','arial','fontsize',12);
set(gca,'fontname','arial','fontsize',12);