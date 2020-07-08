%

clear all
close all
clc

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

%% optimization of log-likelihood
%{
begin_datenum=datenum(begin_date,'yyyy-mm-dd');
end_datenum=datenum(end_date,'yyyy-mm-dd');
ind_begin=find(datetime>=begin_datenum,1,'first');
ind_end=find(datetime<=end_datenum,1,'last');
h=gage_height(ind_begin:ind_end);
discharge=discharge(ind_begin:ind_end);
log_Q_obs=log(discharge);
hmin=min(h); hmax=max(h);
lambda=0.6;           % segmentation rate (number of rating-curve segments per unit length)
nsamp=1;            % number of samples to be drawn
aspace='log';       % space of multiplier parameter ('log' or 'arithmetic')
amin=0.1; amax=5;     % minimum and maximum values of multiplier parameter in appropriate space
bmin=0.5; bmax=3.5; % minimum and maximum values of exponent parameters
alpha=2; beta=0.1;    % parameter of inverse-gamma distribtuion to draw sample from sigma2
h0_min=-5;          % minimum value of first cease-to-flow parameter (in m)
h0_max=hmin;        % maximum value of first cease-to-flow parameter (in m)


nopt = 100;  % number of optimization algorithms
m_max = floor(lambda*(hmax-hmin));

opts = saoptimset('MaxFunEvals',200000,'TolFun',10^-12);
zero_column = zeros(40,1);
count=0;

for m = 1:m_max;
    objfun = @(x)(rc_likeli_opt(x,m,h,log_Q_obs,aspace));
    lb = [h0_min,hmin*ones(1,m-1),1000,h0_min*ones(1,m-1),amin,bmin*ones(1,m),0.00001*ones(1,m)];
    ub = [h0_max,hmax*ones(1,m-1),1000,hmax*ones(1,m-1),amax,bmax*ones(1,m),60*ones(1,m)];
    for opt_ind = 1:nopt
        count = count+1;
        theta0 = joint_priorrnd_opt(m,lambda,hmin,hmax,amin,amax,bmin,bmax,alpha,beta,h0_min,h0_max,aspace,nsamp);
        theta0 = theta0(2:4*m+2);
        theta_opt  = simulannealbnd(objfun,theta0,lb,ub,opts);
        theta(count,:) = [m,add_columns(theta_opt',zero_column)'];
        log_likeli(count)=rc_likeli(theta(count,:),h,log_Q_obs,aspace);
    end
end

% create a probability distribtuion out of optimal theta values using
% log likelihoods as probabilities
theta(log_likeli==inf,:)=[];
log_likeli(log_likeli==inf)=[];
max_log_likeli = max(log_likeli);
diff_log_likeli = log_likeli - max_log_likeli;
probs = exp(diff_log_likeli)/sum(exp(diff_log_likeli));
save('optimal_parameters.mat','theta','probs');
%}
%% draw samples from posterior distribution
%
% open the mat file containing optimal points
filename = fullfile(direc,'optimal_parameters.mat');
load(filename);
samp_space = (1:length(probs))';
probs = probs';

begin_datenum=datenum(begin_date,'yyyy-mm-dd');
end_datenum=datenum(end_date,'yyyy-mm-dd');
ind_begin=find(datetime>=begin_datenum,1,'first');
ind_end=find(datetime<=end_datenum,1,'last');
h=gage_height(ind_begin:ind_end);
discharge=discharge(ind_begin:ind_end);
log_Q_obs=log(discharge);
hmin=min(h); hmax=max(h);
lambda=0.6;           % segmentation rate (number of rating-curve segments per unit length)
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
nsamples=10000;               % number MH samples to be drawn
nchains=8;
thining = 10;

smpl = zeros(nsamples,40,nchains);
theta0 = zeros(nchains,40);
accept = zeros(nchains,1);
for chain_ind=1:nchains
    samp_ind = prob_discrete(samp_space,probs);
    theta0(chain_ind,:) = theta(samp_ind,:);            % seed parameter
    [smpl(:,:,chain_ind),accept(chain_ind)] = mhsample(theta0(chain_ind,:),nsamples,'logpdf',logpdf,'logproppdf',logproppdf,'proprnd',proprnd,'thin',thining);
end


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
uni_smpl=unique(smpl(:,:,1),'rows','stable');
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
        Q_sim_low(ind) = exp(log_Q_sim_low(ind)) - 3*sqrt(sig2); Q_sim_low(Q_sim_low<0) = 0;
        Q_sim_up(ind) = exp(log_Q_sim_up(ind)) + 3*sqrt(sig2);
        
    end
    
    scatter(h,exp(log_Q_sim),'filled'); hold on
    scatter(h,exp(log_Q_obs),'filled','facecolor','r'); hold on
    scatter(h,Q_sim_low,'filled','facecolor','g'); hold on
    scatter(h,Q_sim_up,'filled','facecolor','g'); hold on
%     set(gca,'xscale','log','yscale','log')
    pause(1);
    close all
end
%}

%% R-diagnostic (Gelman and Rubin, 1992)
% compute log_likelihood for each parameter in each chain
%
log_likeli = zeros(nsamples,nchains);
for chain_ind = 1:nchains
    for smp = 1:nsamples
        log_likeli(smp,chain_ind) = rc_likeli(smpl(smp,:,chain_ind),h,log_Q_obs,aspace);
    end
end

% R-diagnostic
sqrt_R = R_diag(log_likeli);
%}