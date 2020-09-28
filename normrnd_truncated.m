% this routine drawn random samples from univarite normal distribution truncated at
% zero
% inputs: mu  = mean
%         sig = standard deviation
%         n   = number of samples to be drawn
% output: r = n samples from truncated normal distribution

function r = normrnd_truncated(mu,sig,n)
    % draw samples from Gaussian distribution
    r = normrnd(mu,sig,[n,1]);
    r = r(r>0);
    
    while length(r)<n
        
        ntmp = n - length(r);
        r = [r;normrnd(mu,sig,[ntmp,1])];
        r = r(r>0);
        
    end
end
