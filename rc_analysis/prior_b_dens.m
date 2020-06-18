% desnity of sample from prior distribution of b
% inputs: bsamp = value of b at which density is ot be computed
%         bmin = minimum value of b
%         bmax = maximum value of b
function dens=prior_b_dens(bsamp,bmin,bmax)
    
    % uniform prior
    n=length(bsamp);
    if sum(bsamp>=bmin)==n && sum(bsamp<=bmax)==n
        dens=prod(1./(bmax-bmin));
    else
        dens=0;
    end
    
end