% error checks

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
amin=0; amax=8;     % minimum and maximum values of multiplier parameter in appropriate space
bmin=0.5; bmax=2.5; % minimum and maximum values of exponent parameters
alpha=4; beta=50;    % parameter of inverse-gamma distribtuion to draw sample from sigma2
h0_min=-5;          % minimum value of first cease-to-flow parameter (in m)
h0_max=hmin;        % maximum value of first cease-to-flow parameter (in m)

%% Imaginery flow values
%{
theta0=joint_priorrnd(lambda,hmin,hmax,amin,amax,bmin,bmax,alpha,beta,h0_min,h0_max,aspace,nsamp);
proprnd=@(x)propMCMCrnd(x,lambda,hmin,hmax,amin,amax,bmin,bmax,alpha,beta,h0_min,h0_max,nsamp);                             % draws samples from proposal distribution 
logproppdf=@(x,y)propMCMCpdf(x,y,lambda,hmin,hmax,amin,amax,bmin,bmax,alpha,beta,h0_min,h0_max);
ind=1;
while ind<inf
    
    theta=propMCMCrnd(theta0,lambda,hmin,hmax,amin,amax,bmin,bmax,alpha,beta,h0_min,h0_max,nsamp);
    log_likeli=rc_likeli(theta,h,log_Q_obs,aspace);
end
%}
%% testing fo faster implementation of rc_est
theta0=joint_priorrnd(lambda,hmin,hmax,amin,amax,bmin,bmax,alpha,beta,h0_min,h0_max,aspace,nsamp);
proprnd=@(x)propMCMCrnd(x,lambda,hmin,hmax,amin,amax,bmin,bmax,alpha,beta,h0_min,h0_max,nsamp);                             % draws samples from proposal distribution 
logproppdf=@(x,y)propMCMCpdf(x,y,lambda,hmin,hmax,amin,amax,bmin,bmax,alpha,beta,h0_min,h0_max);
count=0;
while ind<300
    count=count+1;
    theta=propMCMCrnd(theta0,lambda,hmin,hmax,amin,amax,bmin,bmax,alpha,beta,h0_min,h0_max,nsamp);
    m           =    theta(1);                     % number of segments
    h01         =    theta(2);                     % cease-to-flow parameter of the first segment
    h_s         =    theta(2:m+2);                 % list fo break points (h01 is included in the list of break-points)
    h0_list     =    [h01,theta(m+3:2*m+1)];       % list of cease-to-flow parameters
    a1          =    theta(2*m+2);                 % list of multiplier parameters
    b_list      =    theta(2*m+3:3*m+2);           % list of exponent parameters
    sigma2_list =    theta(3*m+3:4*m+2); 
    
    log_Q1 = rc_est1(h,h_s,a1,b_list,h0_list,aspace);
    log_Q = rc_est(h,h_s,a1,b_list,h0_list,aspace);
    
    scatter(log_Q,log_Q1)
    hold on
    ulimit=max([log_Q;log_Q1]);
    llimit=min([log_Q;log_Q1]);
    plot([llimit ulimit],[llimit ulimit],'color','black')
    pause(0.5);
    close;
end
