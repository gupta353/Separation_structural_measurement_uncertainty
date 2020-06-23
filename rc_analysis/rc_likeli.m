% this routine compute log likelihood of parameters of piecewise power-law
% rating curve
% inputs: theta = parameters of the model
%         h = heights at whoch discharges are to be estimated
%         log_Q_obs = logarithm of observed streamflow
%         aspace = if a is in logarithmic space this parameter should be
%         equal to 'log', otherwise it should be equal to 'arithmetic'
% outputs: log_likeli = log-likelhood of the parameter theta
% Ref: Reitan and Petersen-Overleir (2009). Bayesian methods for estimating
% multi-segment discharge rating curve

function log_likeli=rc_likeli(theta,h,log_Q_obs,aspace)
    
    m           =    theta(1);                     % number of segments
    h01         =    theta(2);                     % cease-to-flow parameter of the first segment
    h_s         =    theta(2:m+2);                 % list fo break points (h01 is included in the list of break-points)
    h0_list     =    [h01,theta(m+3:2*m+1)];       % list of cease-to-flow parameters
    a1          =    theta(2*m+2);                 % list of multiplier parameters
    b_list      =    theta(2*m+3:3*m+2);           % list of exponent parameters
    sigma2_list =    theta(3*m+3:4*m+2);           % list of homoscedastic variance of residuals
    
    log_Q_est = rc_est(h,h_s,a1,b_list,h0_list,aspace);
    res = rc_res(log_Q_obs,log_Q_est);
    
    n = length(log_Q_obs);
    u_brk = h_s(2:end);
    log_likeli=0;
    for seg_ind=1:m
        res_tmp = res(h<u_brk(seg_ind));
        sigma2_tmp = sigma2_list(seg_ind);
        log_likeli = log_likeli - n/2*log(2*pi) - n/2*log(sigma2_tmp) - 1/2/sigma2_tmp*sum(res_tmp.^2);
    end
end