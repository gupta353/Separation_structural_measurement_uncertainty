% this routine plots various uncertainty bands

clear all
close all
clc

direc = 'D:/Research/Thesis_work/Structural_vs_measurement_uncertainty/matlab_codes';

%{
% read rating curve uncertainty data
fname = 'rc_strm_samples_04180500_resampled.mat';
filename = fullfile(direc,'huc_04100003/results',fname);
load(filename)

strm_rc = strm_samps;
clear strm_samps


% read ML uncertainty data
direc_ML = 'E:/covmat_Jefferys_prior';
alpha = 100000;
fname = ['strm_samples_alpha=',num2str(alpha),'.mat'];
filename = fullfile(direc_ML,fname);
load(filename);
strm_ML = strm_samps;
clear strm_samps

% read uncertainty bound obtained by runoff-coefficient method
dist_thresh = 0.20;
fname = ['runoff_coeff_uncer_dist_thresh_',num2str(dist_thresh),'.mat'];
filename = fullfile(direc,'huc_04100003/results/runoff_coeff_uncetainty_bounds_04180500',fname);
load(filename)
% read data in 01-01-2008 to 31-12-2010 range
datenums = datenum('01-01-2005','dd-mm-yyyy'):datenum('31-12-2016','dd-mm-yyyy');
begin_datenum = datenum('01-01-2008','dd-mm-yyyy');
end_datenum = datenum('31-12-2010','dd-mm-yyyy');
ind = find(datenums>=begin_datenum & datenums<=end_datenum);
ub_streamflow_reconstruct = ub_streamflow_reconstruct(ind);
lb_streamflow_reconstruct = lb_streamflow_reconstruct(ind);

% plot data
f1=area(max(strm_ML),'facecolor',[0.2,0.8,0.2],'edgecolor',[0.2,0.9,0.2]);
hold on
f2=area(min(strm_ML),'facecolor',[1,1,1],'edgecolor',[1,1,1]);
f3=plot(min(strm_rc),'--');
f4=plot(max(strm_rc),'--');
f5 = plot(lb_streamflow_reconstruct,'r');
f6 = plot(ub_streamflow_reconstruct,'r');

xlabel('Day since 1 Jan 2008','fontname','arial','fontsize',12);
ylabel('Streamflow (m^3 s^{-1})','fontname','arial','fontsize',12);
set(gca,'xtick',350:10:550,'fontname','arial','fontsize',12,'plotboxaspectratio',[2,1,1])
xlim([400 460])
legend([f1,f3,f5],{'ML','rating curve analysis','runoff-coefficient'},'fontsize',12,'plotboxaspectratio',[2,1,1],'Location','Northwest');
legend('boxoff');

sname = ['uncertainty_ML_rc_runc_dist=',num2str(dist_thresh),'.svg'];
filename = fullfile(direc,'huc_04100003/plots/uncertainty_ML_rc_runc',sname);
fig2svg(filename);
%}

% read ML uncertainty data
%
direc_ML = 'E:/covmat_Jefferys_prior';
datenums = datenum('01-01-2008','dd-mm-yyyy'):datenum('31-12-2008','dd-mm-yyyy');
alphas = [1,50,100000];
colors = {'b','g','r'};
line_pattern = {'-','--','.'};
% plot data
for ind = 1:length(alphas)
    
    fname = ['strm_samples_alpha=',num2str(alphas(ind)),'.mat'];
    filename = fullfile(direc_ML,fname);
    load(filename);
    f(ind)=plot(max(strm_samps),'color',colors{ind},'linewidth',1,'linestyle',line_pattern{ind});
    hold on
    plot(min(strm_samps),'color',colors{ind},'linewidth',1,'linestyle',line_pattern{ind})
    legend_labels{ind} = ['\alpha = ', num2str(alphas(ind))];
    
end

xlabel('Day since 1 Jan 2008','fontname','arial','fontsize',12);
ylabel('Streamflow (m^3 s^{-1})','fontname','arial','fontsize',12);
set(gca,'xtick',350:10:550,'fontname','arial','fontsize',12,'plotboxaspectratio',[2,1,1])
xlim([400 460])
legend(f,legend_labels,'fontname','arial','fontsize',12)
legend('boxoff')
break
sname = 'uncertainty_ML_different_alphas.svg';
filename = fullfile(direc,'huc_04100003/plots',sname);
fig2svg(filename);
%}
