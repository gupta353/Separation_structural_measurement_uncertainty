% draws samples from prior distribution of h0j (cease-to-flow parameters) given the break points
% inputs: h0_min = minimum value of cease-to-flow parameter
%         h_s = break points
%         nsamp = number of samples to be drawn
% outputs: h0j_samp = drawn value of cease-to-flow parameters
%          dens = density of the drawn sample(s)

function [h0j_samp,dens]=prior_h0j(h0_min,h0_max,h_s,nsamp)
    
    % uniform prior
    h_s(1)=[]; h_s(end)=[]; % delete h01 and 'upper limit of the last segment' from h_s
    
    if ~isempty(h_s)
        llimit=h0_min*ones(length(h_s),1);
        ulimit=h_s;
        for samp_ind=1:nsamp
            h0j_samp(:,samp_ind)=unifrnd(llimit,ulimit);
        end
        
        dens=1/prod(ulimit-llimit);
    else
        h0j_samp=[];
        dens=1;
    end
end