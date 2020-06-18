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
    %{
    dens=1;
    %}
    
    % discrete uniform prior
    kmax=floor(lambda*(hmax-hmin));
    if msamp<=kmax && msamp>=1
        dens=1/kmax;
    else
        dens=0;
    end
end