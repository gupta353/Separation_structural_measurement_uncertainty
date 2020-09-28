% this routine samples from a discrete probability distribution
% inputs: samp = elements of sample space (column vector)
%         p    = probabilities of elements contained in sample space (column vector)
% outputs: r  = random variable drawn from the probability space

function r = prob_discrete(samp,p)
    
    if size(samp,2)>1
        error('inputs should be column vectors')
    end
    
    data = [samp,p];
    data = sortrows(data);
    data(:,2) = cumsum(data(:,2));
    
    u = unifrnd(0,1); % uniform random variable
    ind = find(data(:,2)>=u,1,'first');
    r = data(ind,1);
end