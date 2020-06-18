% density of samples drawn from prior_sigma2
% inputs: sigma2_samp = value of drawn sample at which probability density
%                       is to be computed
%         alpha = parameter of (inverse)gamma distribution
%         beta = parameter of (inverse)gamma distribution
% outputs: dens = desnity of samples to be drawn

function dens=prior_sigma2_dens(sigma2_samp,alpha,beta)
    
    % inverse gamma prior
    dens=(beta^alpha*gamma(alpha)*(sigma2_samp).^(alpha+1).*exp(1/beta/sigma2_samp))^-1;
    
end