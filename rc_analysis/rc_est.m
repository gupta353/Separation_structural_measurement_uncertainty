% This routine specifies a multi-segment rating curve where each segment is
% modeled as power-law
% inputs: h=gage heights at which log of streamflow values are to be estimated (in m)
%         a1 =  multiplier paramters for different segments of
%         the rating curve
%         b_list = list of exponent paramters for different segment of
%         the rating curve
%         h0_list = list of cease-to-flow parameters for different segments
%         of the rating curve (in m)
%         h_s = list of break points in gage height (m) (includes h_01 as the first element, a total of n+1 elements where n is the total number of segments)
%         aspace = if parameter 'a' is provided in logarithmic space this,
%         parameter is equal to 'log'; otherwise it should be equal to
%         'arithmetic'
% output: log_Q = logarithm of Q (Q in m^3 s^{-1})

function log_Q = rc_est(h,h_s,a1,b_list,h0_list,aspace)
    
    if strcmp(aspace,'arithmetic')
        a1=log(a1);
    end
    m=length(b_list); % number of rating curve segments
    
    if m>1
        a_list = rc_continuity(a1,h_s,h0_list,b_list);
    else
        a_list=a1;
    end
    % compute log of streamflow
%     u_brk=h_s(2:end);   % upper break point of each rating curve segment
%     for h_ind=1:length(h)
%         
%         htmp=h(h_ind);
%         
%         if htmp>h0_list(1)                  % check if htmp is greater than first cease-to-flow parameter
%             seg=find(htmp<u_brk,1,'first'); % identify the segment
%             log_a=a_list(seg);
%             b=b_list(seg);
%             h0=h0_list(seg);
%             log_Q(h_ind,1)=log_a+b*log(htmp-h0);
%         else
%             log_Q=-inf*ones(length(h),1);
%         end
%     end
    u_brk=h_s(1:end);   % upper break point of each rating curve segment
    log_Q=zeros(length(h),1);
    for seg=2:length(u_brk)
        ind = find(h<u_brk(seg) & h>=u_brk(seg-1));
        htmp = h(ind);
        log_a = a_list(seg-1);
        b = b_list(seg-1);
        h0 = h0_list(seg-1);
        log_Q(ind) = log_a+b*log(htmp-h0);
    end
    log_Q(h<h0_list(1)) = -inf;
    
    if ~isreal(log_Q)
        error('The parameter set generated imaginery flow values. Please check if the parameters are physically consistent')
    end
end