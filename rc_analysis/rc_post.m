% 

% clear all
% close all
% clc

direc='D:/Research/Thesis_work/Structural_vs_measurement_uncertainty/matlab_codes';
begin_date='2001-01-01';    % begin date in yyyy-mm-dd format
end_date='2013-12-31';      % end date in yyyy-mm-dd format

%% rating curvees available on USGS website
%
fname='rc_04180500.txt';
filename=fullfile(direc,'huc_04100003','streamflow_data','rating_curve',fname);
data=readtable(filename,'delimiter','\t');
gage_height=0.3048*data.(1);                  % in m
discharge=0.028316847*data.(2);                    % in cms
date=data.(3);

datetime_wrapper=@(x)datenum(x,'mm/dd/yyyy');
datetime=cellfun(datetime_wrapper,date);

% remove outliers
discharge(gage_height>100)=[];
date(gage_height>100)=[];
datetime(gage_height>100)=[];
gage_height(gage_height>100)=[];

%% draw samples from posterior distribution
%
begin_datenum=datenum(begin_date,'yyyy-mm-dd');
end_datenum=datenum(end_date,'yyyy-mm-dd');
ind_begin=find(datetime>=begin_datenum,1,'first');
ind_end=find(datetime<=end_datenum,1,'last');
h=gage_height(ind_begin:ind_end);
discharge=discharge(ind_begin:ind_end);
log_Q_obs=log(discharge);
hmin=min(h); hmax=max(h);
lambda=0.4;           % segmentation rate (number of rating-curve segments per unit length)
nsamp=1;            % number of samples to be drawn
aspace='log';       % space of multiplier parameter ('log' or 'arithmetic')
amin=0.1; amax=5;     % minimum and maximum values of multiplier parameter in appropriate space
bmin=0.5; bmax=3.5; % minimum and maximum values of exponent parameters
alpha=2; beta=0.1;    % parameter of inverse-gamma distribtuion to draw sample from sigma2
h0_min=-5;          % minimum value of first cease-to-flow parameter (in m)
h0_max=hmin;        % maximum value of first cease-to-flow parameter (in m)


proprnd=@(x)propMCMCrnd(x,lambda,hmin,hmax,amin,amax,bmin,bmax,alpha,beta,h0_min,h0_max,nsamp);                             % draws samples from proposal distribution 
logproppdf=@(x,y)propMCMCpdf(x,y,lambda,hmin,hmax,amin,amax,bmin,bmax,alpha,beta,h0_min,h0_max);                                   % log of transition pdf for proposal distribtuion
logpdf=@(x)(rc_likeli(x,h,log_Q_obs,aspace)+joint_priorpdf(x,lambda,hmin,hmax,amin,amax,bmin,bmax,alpha,beta,h0_min,h0_max));       % lof of target distribution
nsamples=100000;               % number MH samples to be drawn
nchains=8;
thining = 10;

smpl=zeros(nsamples,40,nchains);
for chain_ind=1:nchains
    theta0(:,:,chain_ind) = joint_priorrnd(lambda,hmin,hmax,amin,amax,bmin,bmax,alpha,beta,h0_min,h0_max,aspace,nsamp);            % seed parameter
    [smpl(:,:,chain_ind),accept(chain_ind)] = mhsample(theta0(:,:,chain_ind),nsamples,'logpdf',logpdf,'logproppdf',logproppdf,'proprnd',proprnd,'thin',thining);
end
}

%% compare observed and simulated rating curves
%{
uni_smpl=unique(smpl,'rows','stable');
for theta_ind=1:size(uni_smpl,1)
    theta=uni_smpl(theta_ind,:);
    m       =    theta(1);                     % number of segments
    h01     =    theta(2);                     % cease-to-flow parameter of the first segment
    h_s     =    theta(2:m+2);                 % list fo break points (h01 is included in the list of break-points)
    h0_list =    [h01,theta(m+3:2*m+1)];       % list of cease-to-flow parameters
    a1      =    theta(2*m+2);                 % list of multiplier parameters
    b_list  =    theta(2*m+3:3*m+2);           % list of exponent parameters
    sigma2_list  =    theta(3*m+3:4&m+2);
    
    log_Q_sim = rc_est(h,h_s,a1,b_list,h0_list,aspace);
    scatter(exp(log_Q_obs),exp(log_Q_sim),'filled'); hold on
    llimit = min([exp(log_Q_sim);exp(log_Q_obs)]);
    ulimit = max([exp(log_Q_sim);exp(log_Q_obs)]);
    plot([llimit ulimit],[llimit ulimit],'color','black','linewidth',2)
    pause(0.5);
    close all
end
%}
%% plot the log-log linear relationship for observed and predicted flow values
%{
uni_smpl=unique(smpl,'rows','stable');
for theta_ind=1:size(uni_smpl,1)
    theta=uni_smpl(theta_ind,:);
    m       =    theta(1);                     % number of segments
    h01     =    theta(2);                     % cease-to-flow parameter of the first segment
    h_s     =    theta(2:m+2);                 % list fo break points (h01 is included in the list of break-points)
    h0_list =    [h01,theta(m+3:2*m+1)];       % list of cease-to-flow parameters
    a1      =    theta(2*m+2);                 % list of multiplier parameters
    b_list  =    theta(2*m+3:3*m+2);           % list of exponent parameters
    sigma2_list  =    theta(3*m+3:4*m+2);
    
    log_Q_sim = rc_est(h,h_s,a1,b_list,h0_list,aspace);
    for seg = 1:m
        
        ind=find(h<=h_s(seg+1) & h>=h_s(seg));
        htmp = h(ind);
        h0 = h0_list(seg);      
        
        scatter(htmp-h0,exp(log_Q_sim(ind)),'filled','facecolor','r'); hold on
        scatter(htmp-h0,exp(log_Q_obs(ind)),'filled','facecolor','b'); hold on
    end
    
    for seg = 2:m
        plot([h_s(seg) h_s(seg)],[min(exp(log_Q_sim)) max(exp(log_Q_sim))],'color','black','linewidth',2);
    end
    
    set(gca,'xscale','log','yscale','log')
    pause(0.5);
    close;
end
%}
%% plot the posterior rating curves
%{
uni_smpl=unique(smpl,'rows','stable');
for theta_ind=1:size(uni_smpl,1)
    theta=uni_smpl(theta_ind,:);
    m       =    theta(1);                     % number of segments
    h01     =    theta(2);                     % cease-to-flow parameter of the first segment
    h_s     =    theta(2:m+2);                 % list fo break points (h01 is included in the list of break-points)
    h0_list =    [h01,theta(m+3:2*m+1)];       % list of cease-to-flow parameters
    a1      =    theta(2*m+2);                 % list of multiplier parameters
    b_list  =    theta(2*m+3:3*m+2);           % list of exponent parameters
    sigma2_list  =    theta(3*m+3:4*m+2);
    
    log_Q_sim = rc_est(h,h_s,a1,b_list,h0_list,aspace);
    log_Q_sim_low = log_Q_sim;
    log_Q_sim_up = log_Q_sim;
    for seg = 1:m
        
        ind = find(h<=h_s(seg+1) & h>=h_s(seg));
        sig2 = sigma2_list(seg);
        Q_sim_low(ind) = exp(log_Q_sim_low(ind)) - 3*sqrt(sig2);
        Q_sim_up(ind) = exp(log_Q_sim_up(ind)) + 3*sqrt(sig2);
        
    end
    
    scatter(h,exp(log_Q_sim),'filled'); hold on
    scatter(h,exp(log_Q_obs),'filled','facecolor','r'); hold on
    scatter(h,Q_sim_low,'filled','facecolor','g'); hold on
    scatter(h,Q_sim_up,'filled','facecolor','g'); hold on
    set(gca,'xscale','log','yscale','log')
    pause;
    close all
end
%}

%% R-diagnostic (Gelman and Rubin, 1992)
% compute log_likelihood for each parameter in each chain
for chain_ind = 1:nchains
    for smp = 1:nsamples
        log_likeli(smp,chain_ind) = rc_likeli(smpl(smp,:,chain_ind),h,log_Q_obs,aspace);
    end
end

% compute average at each time-step
count=0;
for t = 20:10:nsamples
    count = count + 1;
    ind = (t/2+1):t;
    average(count,:) = mean(log_likeli(ind,:));         % within-sequence average
    average2(count,:) = average(count,:).^2;            % square of within-sequence average 
    s2(count,:) = var(log_likeli(ind,:));               % within-sequence variance
    n(count,1) = length(ind);                           
    tmp = cov(s2(count,:),average2(count,:));           
    cov1(count,1) = tmp(1,2);
    tmp = cov(s2(count,:),average(count,:));
    cov2(count,1) = tmp(1,2);
    mu_hat(count,1) = mean(average(count,:));
end

W = mean(s2,2);                                 % average of within-sequence variance
var_W  = var(s2,0,2)/nchains;                   % variance of W
B_by_n = var(average,0,2);                        % varaince of within-sequence mean 
var_B  = 2*B_by_n.^2/(nchains-1).*(n.^2);       % variance of B (= B_by_n*n)
cov_W_B = n/nchains.*(cov1 - 2*mu_hat.*cov2);   % covariance of W and B

var_V_hat =  ((n-1)./n).^2.*var_W + ...
    ((nchains+1)./n./nchains).^2.*var_B + ...
    2*(nchains-1)*(n-1)./n.^2/nchains.*cov_W_B;

V_hat = (n-1)./n.*W + (1+1./nchains).*(B_by_n); % Estimate of variance
df = 2*(V_hat.^2)./var_V_hat;                   % Degrees of freedom

R = V_hat./W.*(df./(df-2));                     % R diagnostics
sqrt_R = R.^0.5;                                % Potential scale reduction factor