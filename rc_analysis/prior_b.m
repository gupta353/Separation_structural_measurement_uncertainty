% draws samples from prior distribution of expoent parameter b
% inputs: bmin = minimum value of b
%         bmax = maximum value of b
%         nsamp = number of samples to be drawn
% outputs: bsamp = samples drawn from the distribution of b
%          dens = density of drawn sample

function [bsamp,dens]=prior_b(bmin,bmax,nsamp)
    
    % uniform prior
    bsamp=unifrnd(bmin,bmax,nsamp,1);
    dens=1/(bmax-bmin)*ones(nsamp,1);
    
end