% This routine computes potential scale reduction factor
% inputs: X = A matrix containing multiple sequences of scalar quantity for
%             which Rd has to be computed (each column contains a separate sequence)
%         sqrt_R = Potential scale reduction factor 
% Ref: Gelman and Rubin (1992).

function sqrt_R = R_diag(X)
    
    nsamples = size(X,1);
    nchains = size(X,2);
    count=0;
    for t = 20:10:nsamples
        count = count + 1;
        ind = (t/2+1):t;
        average(count,:) = mean(X(ind,:));         % within-sequence average
        average2(count,:) = average(count,:).^2;            % square of within-sequence average
        s2(count,:) = var(X(ind,:));               % within-sequence variance
        n(count,1) = length(ind);
        tmp = cov(s2(count,:),average2(count,:));
        cov1(count,1) = tmp(1,2);
        tmp = cov(s2(count,:),average(count,:));
        cov2(count,1) = tmp(1,2);
        mu_hat(count,1) = mean(average(count,:));
    end
    
    W = mean(s2,2);                                 % average of within-sequence variance
    var_W  = var(s2,0,2)/nchains;                   % variance of W
    B_by_n = var(average,0,2);                        % varaince of within-sequence mean
    var_B  = 2*B_by_n.^2/(nchains-1).*(n.^2);       % variance of B (= B_by_n*n)
    cov_W_B = n/nchains.*(cov1 - 2*mu_hat.*cov2);   % covariance of W and B
    
    var_V_hat =  ((n-1)./n).^2.*var_W + ...
        ((nchains+1)./n./nchains).^2.*var_B + ...
        2*(nchains-1)*(n-1)./n.^2/nchains.*cov_W_B;
    
    V_hat = (n-1)./n.*W + (1+1./nchains).*(B_by_n); % Estimate of variance
    df = 2*(V_hat.^2)./var_V_hat;                   % Degrees of freedom
    
    R = V_hat./W.*(df./(df-2));                     % R diagnostics
    sqrt_R = R.^0.5;                                % Potential scale reduction factor
    
    plot(20:10:nsamples,sqrt_R);
end