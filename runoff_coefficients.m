% This routine computes runoff coeffcients of a drainage area for all the
% rainfall-runoff events in the user-supplied rainfall-runoff data
% Note: (1) Missing values in rainfall data should be designated as NaNs in
%           input rainfall data
%       (2) Appearance of negative streamflow at the begining of streamflow
%           hydrograph indicates that mrc has over-estimated the streamflows in
%           the previous hydrograph, to correct for over-estimation of
%           streamflow appended streamflows in previous hydrographs is adjusted
%           such that negative streamflows disappear
% Def: Gap period: a time-period with negligible rainfall
clear all
close all
clc

direc='D:/Research/Thesis_work/Structural_vs_measurement_uncertainty/matlab_codes';
darea=2764.0831;               % drainage area in km2
prcp_thresh=1;               % depth of rainfall which is assumed to be negligible for the separation of hydrographs (in mm/day)
length_thresh=7;               % minimum separation between the end and begining of two storms (In other words, minimum length of gap period)
%% load streamflow data, rainfall data, and MRC
strm_fname='streamflow_18.txt';
filename=fullfile(direc,'huc_04100003','streamflow_data',strm_fname);
fid=fopen(filename,'r');
data=textscan(fid,'%s%f','delimiter',char(32),'headerlines',1);
fclose(fid);
strm_date=data{1};
strm_vals=data{2}*0.028316847;      % cfs to cms

prcp_fname='rainfall.txt';
filename=fullfile(direc,'huc_04100003',prcp_fname);
fid=fopen(filename,'r');
data=textscan(fid,'%s%s%f','delimiter',char(32),'headerlines',1);
fclose(fid);
prcp_date=data{1};
prcp_time=data{2};
prcp_vals=data{3};
% replace NaNs by means of previous and next day
indnan=find(isnan(prcp_vals));
for nan_ind=1:length(indnan)
    indtmp=[indnan(nan_ind)-1;indnan(nan_ind)+1];
    prcp_vals(indnan(nan_ind))=mean(prcp_vals(indtmp));
end

if isnan(prcp_vals)
    error('precipitation data contains NaNs')
end

% convert rainfall into daily timescale from hourly timescale
prcp_tmp=reshape(prcp_vals,24,length(prcp_vals)/24);
prcp_vals=sum(prcp_tmp)';


mrc_fname='MRC.txt';
filename=fullfile(direc,'huc_04100003','MRC',mrc_fname);
fid=fopen(filename,'r');
mrc=textscan(fid,'%f%f','delimiter','\t','headerlines',1);
fclose(fid);
mrc=mrc{2};           % in cms

%% Identify periods with negligible rainfall (in other words, gap-periods)
ind_prcp_thresh=find(prcp_vals<prcp_thresh);
differ=ind_prcp_thresh(2:end)-ind_prcp_thresh(1:end-1);
ind_differ=find(differ>1);                                                      % indices of gap time-steps
eind_gap_periods=[ind_prcp_thresh(ind_differ);ind_prcp_thresh(end)];            % end index of each gap-period
bind_gap_periods=[ind_prcp_thresh(1);ind_prcp_thresh(ind_differ(1:end)+1)];     % begin index of each gap-period
gap_period=[bind_gap_periods,eind_gap_periods];
length_gap_periods=(eind_gap_periods-bind_gap_periods)+1;                       % length of gap-periods

% remove gap periods that are smaller in length than length_thresh
ind_gp=find(length_gap_periods<length_thresh);
gap_period(ind_gp,:)=[];
length_gap_periods(ind_gp)=[];
%% Identify rainfall periods using gap periods
pseudo_gap_period=[0,0;gap_period];
rain_period(:,1)=pseudo_gap_period(:,2)+1;                                      % begin index of rain period
rain_period(:,2)=[pseudo_gap_period(2:end,1)-1;length(strm_vals)];              % end index of rain period
if rain_period(end,1)>rain_period(end,2)
    rain_period(end,:)=[];
end

%% Identify streamflow period corresponding to each rainfall period
%
strm_tot=[];
for strm_ind=1:size(rain_period,1)-1
    
    bind=rain_period(strm_ind,1);                     % begin index of curerent rainfall period
    eind=rain_period(strm_ind+1,1)-1;                 % index just before the start of next rainfall period
    period{strm_ind}.rain=prcp_vals(bind:eind);
    period{strm_ind}.streamflow=strm_vals(bind:eind);

%     strm_tot=[strm_tot;period{strm_ind}.streamflow];
%     plot(strm_tot);
%     pause;
end
strm_ind=strm_ind+1;
period{strm_ind}.rain=prcp_vals(rain_period(strm_ind,1):end);
period{strm_ind}.streamflow=strm_vals(rain_period(strm_ind,1):end);
%}

%% Complete the recession curve by appending the data from MRC
%
period=complete_hydrograph(period,mrc);
%}
%{
figure;
for per_ind=1:length(period)
plot(period{per_ind}.streamflow);
pause(1);
end
%}
%% Computation of runoff coefficients
[runoff_coeff,strm_tmp_vol,prcp_tmp_vol]=runoff_coff_comp(period,darea);
figure; hist(runoff_coeff);
figure; plot(prcp_tmp_vol,runoff_coeff,'o')

% save data
for run_ind=1:length(runoff_coeff)
    period{run_ind}.runoff_coefficient=runoff_coeff(run_ind);
end
sname='rainfall_runoff_data.mat';
save_filename=fullfile(direc,'huc_04100003',sname);
save(save_filename,'period')