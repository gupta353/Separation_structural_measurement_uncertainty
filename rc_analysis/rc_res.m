% this routine computes residual between observed and estimated logarithm
% of streamflows
% inputs: log_Q_obs = logarithm of observed streamflows
%         log_Q_est = logarithm of estimated streamflows
% output: res = residuals


function res = rc_res(log_Q_obs,log_Q_est)
    
    if length(log_Q_obs)==length(log_Q_est)
        res=log_Q_obs-log_Q_est;
    else
        error('Observed and estimated dicharge must have same time-length')
    end

end