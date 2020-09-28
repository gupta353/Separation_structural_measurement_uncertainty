% This routine plots observed and ML model predicted streamflow for SJRW
% stations

clear all
close all
clc

direc=['D:/Research/Thesis_work/Structural_vs_measurement_uncertainty'...
    '/matlab_codes/huc_04100003'];

nodeobs='18';

% read observed streamflow data
fname=strcat('streamflow_',nodeobs,'.txt');
filename=fullfile(direc,'streamflow_data',fname);
fid=fopen(filename,'r');
data=textscan(fid,'%s%f','headerlines',1,'delimiter',{'\t',char(32)},...
    'MultipleDelimsAsOne',1);
fclose(fid);
date=data{1};
obs_strm=data{2}*0.028316847;

datenum_wrapper=@(x)datenum(x,'yyyy-mm-dd');
datenums=cellfun(datenum_wrapper,date);

% read streamflow data estimated by ML models
fname=strcat('streamflow_',nodeobs,'.mat');
filename=fullfile(direc,'SJRW_ML_results',fname);
est_strm=load(filename);
est_strm=est_strm.MC1_temp*0.028316847;

mean_est_strm=mean(est_strm,2);
est_strm_025=prctile(est_strm',2.5)';
est_strm_025(est_strm_025<0)=0;
est_strm_975=prctile(est_strm',97.5)';
est_strm_005=prctile(est_strm',0.5)';
est_strm_005(est_strm_005<0)=0;
est_strm_995=prctile(est_strm',99.5)';
est_strm_prctile=[obs_strm,est_strm_025,est_strm_975,est_strm_005,est_strm_995];
est_strm_prctile=sortrows(est_strm_prctile);

% plot scatter
lower_limit=min(est_strm_prctile(:));
upper_limit=max(est_strm_prctile(:));

p1=scatter(obs_strm,mean_est_strm,5,'markerfacecolor',[0.2,0.6,1],'markeredgecolor',[0.2,0.6,1]);
hold on
p2=plot(est_strm_prctile(:,1),est_strm_prctile(:,2),'-.','color',[1,0,0]);
plot(est_strm_prctile(:,1),est_strm_prctile(:,3),'-.','color',[1,0,0]);
p3=plot(est_strm_prctile(:,1),est_strm_prctile(:,4),'--','color',[1,0,0]);
plot(est_strm_prctile(:,1),est_strm_prctile(:,5),'--','color',[1,0,0]);
p4=plot([lower_limit upper_limit],[lower_limit upper_limit],'color','black');

box('on');
box.linewidth=2;
xlabel('Obserevd streamflow (m^3 s^{-1})','fontname','arial','fontsize',12);
ylabel('Estimated streamflow (m^3 s^{-1})','fontname','arial','fontsize',12);
set(gca,'fontname','arial','fontsize',12,'xlim',[lower_limit upper_limit],...
    'ylim',[lower_limit upper_limit],box);
legend([p1,p2,p3,p4],{'Mean','95% prediction','99% prediction interval','1:1 line'},...
    'fontname','arial','fontsize',12,'location','northwest');
legend('boxoff');
title(['Node',char(32),nodeobs],'fontname','arial','fontsize',12);

sname=strcat('ML_observed_vs_predicted_',nodeobs,'.jpeg');
save_filename=fullfile(direc,sname);
% fig2svg(save_filename);
print(save_filename,'-r300','-djpeg');