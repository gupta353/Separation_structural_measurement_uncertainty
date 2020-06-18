% proposal distribution for MCMC sampling

clear all
close all
clc


hmin=0.017; hmax=7;
lambda=1;
nsamp=1;
amin=0; amax=8;
bmin=0.5; bmax=2.5;
alpha=3; beta=1;
aspace='log';
h0_min=-5;     % in m
h0_max=hmin;    % in m

theta_old=[5,-4.32,0.95,1.99,4.04,5.66,1000,-1.75,1.70,3.72,-3.32,7.77,...
    2.41,1.47,2.10,0.78,1.34,0.17,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

m_old       =    theta_old(1);                     % number of segments
h01_old     =    theta_old(2);                     % cease-to-flow parameter of the first segment
h_s_old     =    theta_old(2:m_old+2);                 % list fo break points (h01 is included in the list of break-points)
h0_list_old =    [h01_old,theta_old(m_old+3:2*m_old+1)];       % list of cease-to-flow parameters
a1_old      =    theta_old(2*m_old+2);           % list of multiplier parameters
b_list_old  =    theta_old(2*m_old+3:3*m_old+2);           % list of exponent parameters
sigma2_old  =    theta_old(3*m_old+3);

%% draw m (number of segments)
% probabilities (p(m_old)=0.8, p(m_old+1)=0.1,p(m_old-1)0.1)
m_max=floor(lambda*(hmax-hmin));
samp_space=[m_old*ones(8,1);m_old-1;m_old+1]; % population from which a sample can be drawn uniformly
m_new=datasample(samp_space,1);
% make sure that jump is in the support
if m_new>m_max
    m_new=m_old-1;
elseif m_new<1 && m_max>=2
    m_new=m_old+1;
end

%% draw h01
if m_new==m_old
    half_range=min(abs(h01_old-h0_max),abs(h01_old-h0_min));
    sigma_h01=half_range/3;
    h01_new=h01_old+normrnd(0,sigma_h01);
else
    h01_new=prior_h01(h0_min,h0_max,nsamp);
end

%% draw h_s (Gaussian distribution centered around the current value)
if m_new==m_old
    h_s_new=h_s_old; h_s_new(1)=h01_new;  % check
    for hs_ind=2:length(h_s_new)-1
        hs_tmp=h_s_new(hs_ind);
        hs_min=h_s_new(hs_ind-1);
        hs_max=h_s_new(hs_ind+1);
        half_range=min(abs(hs_tmp-hs_min),abs(hs_tmp-hs_max));
        sigma_hs=half_range/3.1;
        
        h_s_new(hs_ind)=hs_tmp+normrnd(0,sigma_hs);
    end
else
    h_s_new=prior_h_s(m_new,hmin,hmax,nsamp)';
end

%% draw h0j
h0j_new=prior_h0j(h0_min,h0_max,h_s_new',nsamp)';
h0j_new=[h01_new,h0j_new];

%% draw exponent parameters
if m_new==m_old
    b_list_new=b_list_old;
    for b_ind=1:length(b_list_new)
        btmp=b_list_old(b_ind);
        half_range=min(abs(btmp-bmin),abs(btmp-bmax));
        sigma_b=half_range/3;
        
        b_list_new(b_ind)=btmp+normrnd(0,sigma_b);
    end
else
    b_list_new=prior_b(bmin,bmax,m_new)';
end

%% draw multiplier parameters
if m_new==m_old
    half_range=min(abs(a1_old-amin),abs(a1_old-amax));
    sigma_a=half_range/3;
    a1_new=a1_old+normrnd(0,sigma_a);
else
    a1_new=prior_a(amin,amax,nsamp);
end

%% draw sigma2
sigma2_new=prior_sigma2(alpha,beta,nsamp);

%% create theta_new
theta_new=[m_new,h01_new,h_s_new(2:end),h0j_new(2:end),a1_new,b_list_new,sigma2_new];
theta_zeros=zeros(1,length(theta_old));
theta_new=add_columns(theta_new',theta_zeros')';