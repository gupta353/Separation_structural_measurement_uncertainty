% this routine computes the log of multiplier parameter (log(a)) of all the
% rating curve segments by the way of continuity principle
% inputs: a1 = multiplier parameter of the first rating-curve segment
%         hs = break point of rating curve segments
%         h0_list = list of cease of flow parameters
%         b_list = list of exponent parameters
%         aspace = if parameter 'a' is provided in logarithmic space this,
%         parameter is equal to 'log'; otherwise it should be equal to
%         'arithmetic'

function a_list = rc_continuity(a1,hs,h0_list,b_list,aspace)
    
    if strcmp(aspace,'arithmetic')
        a1=log(a1);
    end
   
    a_list=[a1;zeros(length(b_list)-1,1)];
    hs(1)=[]; hs(end)=[];
    
    % estimate log of multiplier parameters in iterations
    for a_ind=2:length(a_list)
        
        a_old=a_list(a_ind-1);
        b_old=b_list(a_ind-1);
        b_new=b_list(a_ind);
        h0_old=h0_list(a_ind-1);
        h0_new=h0_list(a_ind);
        hs_iter=hs(a_ind-1);
        
        a_list(a_ind)=a_old + b_old*log(hs_iter-h0_old) - b_new*log(hs_iter-h0_new);
        
    end
    
     if strcmp(aspace,'arithmetic')
        a_list=exp(a_list);
     end
    
end