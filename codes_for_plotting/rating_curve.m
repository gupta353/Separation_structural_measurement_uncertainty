% This routine plots rating curve for a given station in a given
% time-period

clear all
close all
clc

direc='D:/Research/Thesis_work/Structural_vs_measurement_uncertainty/matlab_codes';
begin_date='1984-01-01';    % begin date in yyyy-mm-dd format
end_date='2016-12-31';      % end date in yyyy-mm-dd format

%% load data
fname='rc_04180500.txt';
filename=fullfile(direc,'huc_04100003','streamflow_data',fname);
data=readtable(filename,'delimiter','\t');
gage_height=0.3048*data.(1);                  % in ft
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
begin_datenum=datenum(begin_date,'yyyy-mm-dd');
end_datenum=datenum(end_date,'yyyy-mm-dd');

ind_begin=find(datetime>=begin_datenum,1,'first');
ind_end=find(datetime<=end_datenum,1,'last');
scatter(gage_height(ind_begin:ind_end),discharge(ind_begin:ind_end),'filled');
set(gca,'xscale','log','yscale','log')
hold on

% best-fit line
%
x=log(gage_height(ind_begin:ind_end));
y=log(discharge(ind_begin:ind_end));
% b=exp(x(y==min(y))); b=b(1);
% fit data in a for loop
b=0.1;
rmse_diff_frac=inf;
rmse=0;
while rmse_diff_frac>0.01
    
    x=log(exp(x)-b);
    mdl=fitlm(x,y);
    ypred=mdl.Fitted;
    coeffs=mdl.Coefficients.Estimate;
    b=mean(exp(x)-exp((y-coeffs(1))/coeffs(2)));
    rmse_diff_frac=abs(mdl.RMSE-rmse)/rmse;
    rmse=mdl.RMSE;

end


plot_data=[exp(x),exp(ypred)];
plot_data=sortrows(plot_data);
plot(plot_data(:,1),plot_data(:,2),'r','linewidth',2)
xlabel('gage height (m)','fontname','arial','fontsize',12);
ylabel('Streamflow (m^3 s^{-1})','fontname','arial','fontsize',12);
set(gca,'fontname','arial','fontsize',12)
%}