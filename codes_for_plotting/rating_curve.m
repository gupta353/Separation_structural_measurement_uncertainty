% This routine plots rating curve for a given station in a given
% time-period

clear all
close all
clc

direc='D:/Research/Thesis_work/Structural_vs_measurement_uncertainty/matlab_codes';
begin_date='2009-03-13';    % begin date in yyyy-mm-dd format
end_date='2011-03-05';      % end date in yyyy-mm-dd format

%% rating curvees available on USGS website
%
fname='rc_04180500.txt';
filename=fullfile(direc,'huc_04100003','streamflow_data','rating_curve',fname);
data=readtable(filename,'delimiter','\t');
gage_height=0.3048*data.(1);                  % in m
discharge=0.028316847*data.(2);                    % in cfs
date=data.(3);

datetime_wrapper=@(x)datenum(x,'mm/dd/yyyy');
datetime=cellfun(datetime_wrapper,date);

%% find missing value
ind=find(gage_height>200);
gage_height(ind)=[];
discharge(ind)=[];
datetime(ind)=[];
%% plot rating curve
%{
begin_datenum=datenum(begin_date,'yyyy-mm-dd');
end_datenum=datenum(end_date,'yyyy-mm-dd');

ind_begin=find(datetime>=begin_datenum,1,'first');
ind_end=find(datetime<=end_datenum,1,'last');
scatter(gage_height(ind_begin:ind_end),discharge(ind_begin:ind_end),'filled');
% set(gca,'xscale','log','yscale','log')
grid on
%}

%% rating curve obtaoned through communication with USGS
%
fname_list={'7.1.txt','8.0.txt','8.1.txt','9.0.txt','10.0.txt'};
% figure;
hold on
for fname_ind=1:length(fname_list)
    
    fname=fname_list{fname_ind};
    filename=fullfile(direc,'huc_04100003','streamflow_data','rating_curve',fname);
    fid=fopen(filename,'r');
    data=textscan(fid,'%f%f');
    gage_height=data{2}*0.3048;         % in m
    discharge=data{1}*0.028316847;      % in cms
    
    scatter(gage_height,discharge,'filled');
    
end
legend(fname_list,'fontname','arial','fontsize',12,'location','northwest');
legend('boxoff')
set(gca,'fontname','arial','fontsize',12)
set(gca,'xscale','log','yscale','log')
grid on
xlabel('Gage height (m)','fontname','arial','fontsize',12);
ylabel('Streamflow (m^3 s^{-1})','fontname','arial','fontsize',12);
%}