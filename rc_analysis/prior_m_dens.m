% probability density drawn of a sample drawm from prior_m
% inputs: msamp = number of rating curve segments
%         lambda = average number of segments per-unit length in
%                  gage-height space
%         hmin = minimum recorded gage-height
%         hmax = maximum recorded gage-height
% outputs: dens = probability density of msamp

function dens=prior_m_dens(msamp,lambda,hmin,hmax)
    
    % poisson prior
    %{
    len=hmax-hmin;                    % total length
    lambda_len=lambda*len;
    dens=(lambda_len.^msamp).*exp(-lambda_len)./factorial(msamp)/(1-exp(-lambda_len));
    %}
    
    % fixed value
    %
    if msamp==2
        dens=1;
    else
        dens=0;
    end
    %}
    
    % discrete uniform prior
%     kmax=floor(lambda*(hmax-hmin));
%     if msamp<=kmax && msamp>=1
%         dens=1/kmax;
%     else
%         dens=0;
%     end
    
    % non-uniform prior (p_(i+1)=2*p_i)
%     kmax=floor(lambda*(hmax-hmin));
%     mult=10^3;      % GP multiplier
%      if msamp<=kmax && msamp>=1
%         dens=mult^(msamp-1)/(mult^kmax-1);
%     else
%         dens=0;
%     end
end