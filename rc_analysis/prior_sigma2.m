% drawn samples from prior distribution of \sigma^2
% inputs: alpha = parameter of (inverse)gamma distribution
%         beta = parameter of (inverse)gamma distribution
%         nsamp = number of samples to be drawn
% outputs: sigma2_samp = samples of drawn sigma2
%          dens = desnity of samples to be drawn

function [sigma2_samp,dens]=prior_sigma2(alpha,beta,nsamp)
    
    % inverse gamma prior
    sigma2_inv_samp=gamrnd(alpha,beta,nsamp,1);
    sigma2_samp=1./sigma2_inv_samp;

    dens=beta^alpha/gamma(alpha)*(1./sigma2_samp).^(alpha+1).*exp(-beta./sigma2_samp);
    
end