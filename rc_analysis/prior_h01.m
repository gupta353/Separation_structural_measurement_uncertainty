% draws samples from prior distribution of h01
% inputs: h0_min = minimum value of cease-to-flow parameters
%         h0_max = maximum value of cease-to-flow parameters
%         nsamp = number of samples to be drawn
% outputs: h01_samp = value of drawn h01
%          dens = probability density value of drawn h01

function [h01_samp,dens]=prior_h01(h0_min,h0_max,nsamp)
    
    % uniform prior
    h01_samp=unifrnd(h0_min,h0_max,nsamp,1);
    dens=1/(h0_max-h0_min)*ones(nsamp,1);
    
end