% this routine compute log likelihood of parameters of piecewise power-law
% rating curve
% inputs: theta = parameters of the model
%         m = number of segments
%         h = heights at whoch discharges are to be estimated
%         log_Q_obs = logarithm of observed streamflow
%         aspace = if a is in logarithmic space this parameter should be
%         equal to 'log', otherwise it should be equal to 'arithmetic'
% outputs: log_likeli = log-likelhood of the parameter theta
% Ref: Reitan and Petersen-Overleir (2009). Bayesian methods for estimating
% multi-segment discharge rating curve

function log_likeli=rc_likeli_opt(theta,m,h,log_Q_obs,aspace)
    
    h01         =    theta(1);                     % cease-to-flow parameter of the first segment
    h_s         =    theta(1:m+1);                 % list fo break points (h01 is included in the list of break-points)
    h0_list     =    [h01,theta(m+2:2*m)];       % list of cease-to-flow parameters
    a1          =    theta(2*m+1);                 % list of multiplier parameters
    b_list      =    theta(2*m+2:3*m+1);           % list of exponent parameters
    sigma2_list =    theta(3*m+2:4*m+1);           % list of homoscedastic variance of residuals
    
    log_Q_est = rc_est_opt(h,h_s,a1,b_list,h0_list,aspace);
        
    u_brk = h_s(2:end);
    
    % multiplicative error model
    %{
    res = rc_res(log_Q_obs,log_Q_est);
    log_likeli=0;
    for seg_ind=1:m
        ind = find(h<u_brk(seg_ind));
        res_tmp = res(ind);
        ntmp = length(ind);
        sigma2_tmp = sigma2_list(seg_ind);
        log_likeli = log_likeli - ntmp/2*log(2*pi) - ntmp/2*log(sigma2_tmp) - 1/2/sigma2_tmp*sum(res_tmp.^2);
    end
    %}
    
    % Addtive error model with estimated response truncated at zero
    %
    Q_est = exp(log_Q_est);
    Q_obs = exp(log_Q_obs);
    res = rc_res(Q_obs,Q_est);
    log_likeli=0;
    for seg_ind=1:m
        ind = find(h<u_brk(seg_ind));
        res_tmp = res(ind);
        ntmp = length(ind);
        Q_esttmp = Q_est(ind);
                
        sigma2_tmp = sigma2_list(seg_ind);
        log_likeli = log_likeli - ntmp/2*log(2*pi) - ntmp/2*log(sigma2_tmp) - 1/2/sigma2_tmp*sum(res_tmp.^2) + ...
            log(prod(Indicator(res_tmp,-Q_esttmp,inf(ntmp,1)))) - ...
            log(1-prod(normcdf(0,Q_esttmp,sqrt(sigma2_tmp))));
   end
    %}
    log_likeli = -log_likeli;
end