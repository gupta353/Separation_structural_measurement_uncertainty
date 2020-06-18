% draws samples from prior distribution of h_s conditioned on number of
% segments
% inputs: m = number of rating-curve segments
%         hmin = minimum value of breaking points
%         h_max = maximum value of breaking points
%         nsamp = number of samples to be drawn
% output: h_s_samp = drawn sample(s) of break points (vector or matrix of
%                    column vectors)
%         dens = probability density of drawn sample(s)

function [h_s_samp,dens]=prior_h_s(m,hmin,hmax,nsamp)
    
    % conditional ordered-uniform prior
    if m>1
        
        h_s_samp=unifrnd(hmin,hmax,m-1,nsamp);
        h_s_samp=[hmin*ones(1,nsamp);h_s_samp];
        h_s_samp=sort(h_s_samp,1);
        h_s_samp=[h_s_samp;1000*ones(1,nsamp)];
        
        dens=factorial(m-1)/(hmax-hmin)^(m-1)*ones(1,nsamp);
    else
        h_s_samp=[hmin;1000];
        dens=1;
    end
end