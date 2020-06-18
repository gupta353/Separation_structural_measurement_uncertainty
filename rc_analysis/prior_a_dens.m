% desnity of sample from prior distribution of a
% inputs: asamp = sample of a at which density os to be computed
%         amin = minimum value of (log of) a
%         amax = maximum value of (log of) a
% outputs:dens = density of the drawn sample
% Note: this function is valid only for scalar values of asamp

function dens=prior_a_dens(asamp,amin,amax)
    
    % uniform prior
    a=asamp(1);
    if a>=amin && a<=amax
        dens=1/(amax-amin);
    else
        dens=0;
    end
    
end