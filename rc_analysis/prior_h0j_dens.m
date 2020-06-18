% density of sample drawn from prior_h0j
% inputs: h_0j = value of the parameter at which density is to be computed
%         h0_min = minimum value of cease-to-flow parameter
%         h_s = break points
% outputs: dens = density of the drawn sample(s)

function dens=prior_h0j_dens(h_0j,h0_min,h0_max,h_s)
    
    % uniform prior
    h_s(1)=[]; h_s(end)=[]; % delete h01 and 1000 from h_s
    
    if ~isempty(h_s)
        llimit=h0_min*ones(1,length(h_s));
        ulimit=h_s;
        n=length(h_s);
        if sum(h_0j>=llimit)==n && sum(h_0j<=ulimit)==n
            dens=1/prod(ulimit-llimit);
        else
            dens=0;
        end
    else
        dens=1;
    end
    
end