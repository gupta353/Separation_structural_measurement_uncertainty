% This routine plots histograms of variance and correlation for different
% values of coarsening parameter (\alpha)

clear all
close all
clc

direc='D:/Research/Thesis_work/Structural_vs_measurement_uncertainty/matlab_codes';
alpha_vals=[10,100,500,1000,2000,5000,10000,20000];
colorpal=[0,0,1;0,1,0;1,0,0;1,1,0;1,0,1;0,1,1;0,0,0;0.50,0.5,0];

%% plot variance histograms in a for loop
%{
figure; hold on
for fig_ind=1:length(alpha_vals)
    
    alpha_tmp=alpha_vals(fig_ind);
    fname=strcat('covmat_alpha =',num2str(alpha_tmp),'.mat');
    filename=fullfile(direc,'huc_04100003','results',...
        'covmat_Jefferys_prior',fname);
    load(filename);
    legend_lables{fig_ind}=strcat('\alpha = ',num2str(alpha_vals(fig_ind)));
    
    [f,x]=ksdensity(squeeze(C(1,1,:)));
    plot(x,f,'color',colorpal(fig_ind,:),'linewidth',2)
    clear C
end
xlabel('variance at time-step 1 ((m^3 s^{-1})^2)','fontname','arial','fontsize',12);
ylabel('probability density (1/(m^3 s^{-1})^2))','fontname','arial','fontsize',12);
set(gca,'fontname','arial','fontsize',12);
legend(legend_lables,'fontname','arial','fontsize',12);
legend('boxoff');

% save figure
sname='varaince_t=1.fig';
save_filename=fullfile(direc,'huc_04100003','results',...
        'covmat_Jefferys_prior',sname);
save(save_filename)
%}
%% plot correlation histograms in a for loop
figure; hold on
for fig_ind=1:length(alpha_vals)
    
    alpha_tmp=alpha_vals(fig_ind);
    fname=strcat('covmat_alpha =',num2str(alpha_tmp),'.mat');
    filename=fullfile(direc,'huc_04100003','results',...
        'covmat_Jefferys_prior',fname);
    load(filename);
    legend_lables{fig_ind}=strcat('\alpha = ',num2str(alpha_vals(fig_ind)));
    
    data=squeeze(C(2,1,:))./sqrt(squeeze(C(2,2,:)))./sqrt(squeeze(C(1,1,:)));
    [f,x]=ksdensity(data);
    plot(x,f,'color',colorpal(fig_ind,:),'linewidth',2)
    clear C data
end
xlabel('variance at time-step 1 ((m^3 s^{-1})^2)','fontname','arial','fontsize',12);
ylabel('probability density (1/(m^3 s^{-1})^2))','fontname','arial','fontsize',12);
set(gca,'fontname','arial','fontsize',12);
legend(legend_lables,'fontname','arial','fontsize',12);
legend('boxoff');

% save figure
sname='corr_t=1_2.fig';
save_filename=fullfile(direc,'huc_04100003','results',...
        'covmat_Jefferys_prior',sname);
save(save_filename)