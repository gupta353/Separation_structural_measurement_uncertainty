% This routine plots rating curve for a given station in a given
% time-period

clear all
close all
clc

direc='D:/Research/Thesis_work/Structural_vs_measurement_uncertainty/matlab_codes';
begin_date='2008-01-01';    % begin date in yyyy-mm-dd format
end_date='2013-12-31';      % end date in yyyy-mm-dd format

%% load data
fname='rc_official_filtered_04180500.txt';
filename=fullfile(direc,'huc_04100003','streamflow_data',fname);
data=readtable(filename,'delimiter','\t');
gage_height=0.3048*data.(3);                  % in m
discharge=0.028316847*data.(2);                    % in cfs
date=data.(1);

datetime_wrapper=@(x)datenum(x,'mm/dd/yyyy');
datetime=cellfun(datetime_wrapper,date);

%% find missing value
ind=find(gage_height>200);
gage_height(ind)=[];
discharge(ind)=[];
datetime(ind)=[];
%% plot rating curve
begin_datenum=datenum(begin_date,'yyyy-mm-dd');
end_datenum=datenum(end_date,'yyyy-mm-dd');

ind_begin=find(datetime>=begin_datenum,1,'first');
ind_end=find(datetime<=end_datenum,1,'last');
scatter(gage_height(ind_begin:ind_end),discharge(ind_begin:ind_end),'filled');
set(gca,'xscale','log','yscale','log')
grid on
hold on