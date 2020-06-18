% this function draws sample from joint prior distribution
% inputs: lambda=segmentation rate (averge number of segments per unit length)
%         hmin=minimum gage height in the data
%         hmax=maximum gage height in the data
%         amin=minimum value of the multiplier parameter of the first
%              segment
%         amx=maximum value of the multiplier parameter of the second
%             segment
%         bmin=minimum value of the exponent parameter
%         bmax=maximum value of the exponent parameter
%         alpha and beta=parameter of the inverse gamma distribtuion to be used to
%         draw samples from population of sigma2
%         aspace = 'log' is the multiplier parameter is in log-space,
%                   'arithmetic' otherwise 
%         nsamp = number of samples to be drawn
% outputs: theta = parameters from the joint distribution

function theta=joint_priorrnd(lambda,hmin,hmax,amin,amax,bmin,bmax,alpha,beta,h0_min,h0_max,aspace,nsamp)
        
    msamp = prior_m(lambda,hmin,hmax,nsamp);
    h01_samp = prior_h01(h0_min,h0_max,nsamp);
    hs_samp = prior_h_s(msamp,hmin,hmax,nsamp);
    h0j_samp = prior_h0j(h0_min,h0_max,hs_samp,nsamp);
    bsamp = prior_b(bmin,bmax,msamp);
    asamp = prior_a(amin,amax,nsamp);
    sigma2_samp =prior_sigma2(alpha,beta,nsamp);
    
    theta=[msamp;h01_samp;hs_samp(2:end);h0j_samp;asamp;bsamp;sigma2_samp]'; % first element of hs_samp contains h01_samp
    
    % apppend theta by zeros so that the number of elements remain constant
    theta_zeros=zeros(1,40);
    theta=add_columns(theta_zeros',theta')';
    
%     joint_dens=mdens*h01_dens*hs_dens*h0_dens*prod(adens)*prod(bdens)*sigma2_dens;
    
end