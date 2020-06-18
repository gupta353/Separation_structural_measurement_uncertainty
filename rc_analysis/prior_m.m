% drawn samples from prior distribution of m
% inputs: lambda = average number of segments per-unit length in
%                  gage-height space
%         hmin = minimum recorded gage-height
%         hmax = maximum recorded gage-height
%         nsamp = number of samples to be drawn
% outputs: msamp = number of segments
%          dens = probability density of msamp

function [msamp,dens]=prior_m(lambda,hmin,hmax,nsamp)
    
    % poisson prior
    %{
    len=hmax-hmin;                    % total length
    lambda_len=lambda*len;
    msamp=0;
    while msamp==0
        msamp=poissrnd(lambda_len,nsamp,1);
    end
    dens=(lambda_len.^msamp).*exp(-lambda_len)./factorial(msamp)/(1-exp(-lambda_len));
    %}
    
    % fixed value
    %{
    msamp=3;
    dens=1;
    %}
    
    % discrete uniform prior
    kmax=floor(lambda*(hmax-hmin));
    msamp=randsample(kmax,nsamp);
    dens=1/kmax;
end