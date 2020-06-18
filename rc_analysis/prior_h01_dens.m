% density of samples drawn from prior_h01
% inputs: h01 = value of h01 (cease-to-flow parameter) at which density
%               is to be computed
%         h0_min = minimum value of cease-to-flow parameters
%         h0_max = maximum value of cease-to-flow parameters
% outputs: dens = probability density value of drawn h01
% Note: this funciton is valid only for scalar values of h01

function dens=prior_h01_dens(h01,h0_min,h0_max)
    
    % uniform prior
    if h01>=h0_min && h01<=h0_max
        dens=1/(h0_max-h0_min);
    else
        dens=0;
    end
    
end