% density of sample drawn from prior_h_s
% inputs: h_s = breaking point values
%         m = number of rating-curve segments
%         hmin = minimum value of breaking points
%         hmax = maximum value of breaking points
% output: dens = probability density of drawn sample(s)

function dens=prior_h_s_dens(h_s,m,hmin,hmax)
    
    % conditional uniform prior ordered
    h_s(1)=[];h_s(end)=[];
    n=length(h_s);
    if m>1
        if sum(h_s>=hmin)==n && sum(h_s<=hmax)==n
            dens=factorial(m-1)/(hmax-hmin)^(m-1);
        else
            dens=0;
        end
    else
        dens=1;
    end
    
end