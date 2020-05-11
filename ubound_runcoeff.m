% This routine computes uncertainty bounds over streamflow for each
% rainfall-runoff period
% Ref: Beven (2019). Towards hypothesis testing in inexact sciences

clear all
close all
clc

direc='D:/Research/Thesis_work/Structural_vs_measurement_uncertainty/matlab_codes';

%% load data
fname='rainfall_runoff_data.mat';
filename=fullfile(direc,'huc_04100003',fname);
load(filename);

%% record all runoff-coefficients
for per_ind=1:length(period)
    
    runoff_coeff(per_ind)=period{per_ind}.runoff_coefficient;
    
end
runoff_coeff(runoff_coeff>3)=[];
runoff_coeff_max=max(runoff_coeff);
runoff_coeff_min=min(runoff_coeff);
%% compute uncertainty bounds by assigning equal membership to all runoff coefficients
for per_ind=1:length(period)
    
    period{per_ind}.ub_streamflow=period{per_ind}.completed_streamflow*runoff_coeff_max/period{per_ind}.runoff_coefficient;
    period{per_ind}.lb_streamflow=period{per_ind}.completed_streamflow*runoff_coeff_min/period{per_ind}.runoff_coefficient;
    
end

% plot data
for per_ind=1:length(period)
    plot(period{per_ind}.lb_streamflow);
    hold on; 
    plot(period{per_ind}.ub_streamflow);
    plot(period{per_ind}.completed_streamflow,'r');
    pause;
    hold off;
end

%% compute membership of each runoff coefficient
