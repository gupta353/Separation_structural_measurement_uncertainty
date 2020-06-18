% draws sample(s) from prior distribution of (log) of multiplier parameter a
% of the first segment of rating curve
% inputs: amin = minimum value of (log of) a
%         amax = maximum value of (log of) a
%         nsamp = number of samples to be drawn
% outputs: asamp = drawn sample(s) of a
%          dens = density of the drawn sample(s)

function [asamp,dens]=prior_a(amin,amax,nsamp)
    
    % uniform prior
    asamp=unifrnd(amin,amax,nsamp,1);
    dens=1/(amax-amin)*ones(nsamp,1);
    
end