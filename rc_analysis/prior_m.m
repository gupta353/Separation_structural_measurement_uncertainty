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
%     kmax=floor(lambda*(hmax-hmin));
%     msamp=randsample(kmax,nsamp);
%     dens=1/kmax;
    
    % non-uniform prior (p_(i+1)=2*p_i)
    kmax=floor(lambda*(hmax-hmin));
    mult=10^3;      % GP multilier
    samp_space=[];
    for ind=1:kmax
        samp_space=[samp_space;ind*ones(mult^(ind-1),1)];
    end
    msamp = datasample(samp_space,1);
    
    dens=mult^(msamp-1)/(2^kmax-1);
    
end