% This function returns the log of joint density of the prior distribtuion
% at the sample value drawn from join_priorend.m
% inputs: theta = parameter set from the joint distribution at which
%                 density is to be computed
%         lambda = segmentation rate
%         hmin = minimum gage height in the data
%         hmax = maximum gage height in the data
%         amin = minimum value of the multiplier parameter
%         amax = maximum value of the multiplier parameter
%         bmin = minimum value of the exponent parameter
%         bmax = maximum value of the exponent parameter
%         h0_min and h0_max = minimum and maximum value of cease-to-flow
%                             parameter
%         alpha and beta = a parameter of the inverse gamma distribtuion to be used to
%                          draw samples from population of sigma2
%         nsamp = number of samples to be drawn
% outputs: log_joint_dens = joint density

function log_joint_dens=joint_priorpdf(theta,lambda,hmin,hmax,amin,amax,bmin,bmax,alpha,beta,h0_min,h0_max)
        
        m       =    theta(1);                     % number of segments
        h01     =    theta(2);                     % cease-to-flow parameter of the first segment
        h_s     =    theta(2:m+2);                 % list fo break points (h01 is included in the list of break-points)
        h0_list =    [h01,theta(m+3:2*m+1)];       % list of cease-to-flow parameters
        a1      =    theta(2*m+2);                 % list of multiplier parameters
        b_list  =    theta(2*m+3:3*m+2);           % list of exponent parameters
        sigma2  =    theta(3*m+3);                 % homoscedastic variance of residuals
        
        densm=prior_m_dens(m,lambda,hmin,hmax);
        dens_h01=prior_h01_dens(h01,h0_min,h0_max);
        dens_h_s=prior_h_s_dens(h_s,m,hmin,hmax);
        dens_h_0j=prior_h0j_dens(h0_list(2:end),h0_min,h0_max,h_s);
        densa=prior_a_dens(a1,amin,amax);
        densb=prior_b_dens(b_list,bmin,bmax);
        dens_sigma2=prior_sigma2_dens(sigma2,alpha,beta);
        
        log_joint_dens=log(densm)+log(dens_h01)+log(dens_h_s)+log(dens_h_0j)+log(densa)+log(densb)+log(dens_sigma2);
               
end