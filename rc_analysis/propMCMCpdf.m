% This routine computes density proposed parameter
% Inputs: theta_old = old set of parameters
%         theta_new = propsed set of parameters
%         Rest of the parameters are the same as in propMCMCrnd
% Outputs: log_dens_prop_prior = logarithm of density of the proposed
%                                set of parameters
%% Test example
% hmin=0.3871; hmax=6.1570;
% lambda=1;           % segmentation rate (number of rating-curve segments per unit length)
% aspace='log';       % space of multiplier parameter ('log' or 'arithmetic')
% amin=0; amax=8;     % minimum and maximum values of multiplier parameter in appropriate space
% bmin=0.5; bmax=2.5; % minimum and maximum values of exponent parameters
% alpha=3; beta=1;    % parameter of inverse-gamma distribtuion to draw sample from sigma2
% h0_min=-5;          % minimum value of first cease-to-flow parameter (in m)
% h0_max=hmin;        % maximum value of first cease-to-flow parameter (in m)

% theta_old=[5,-0.35,1.62,5.37,5.96,6.10,1000,-3.20,5.34,1.23,-3.53,1.65,2.23,...
%     1.44,1.52,1.51,2.19,0.34,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% theta_new=[5,-0.48,1.65,5.21,6.06,6.12,1000,-1.20,-2.29,3.37,5.50,1.69,2.15,1.50,1.36,1.40,2.18,0.62,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% theta_new=[4,-2.10,3.41,4.47,5.99,1000,-3.32,-0.60,-1.58,5.54,2.23,0.96,1.79,0.56,0.35,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

function log_dens_prop_prior = propMCMCpdf(theta_old,theta_new,lambda,hmin,hmax,amin,amax,bmin,bmax,alpha,beta,h0_min,h0_max)
    %% cirrent parameters
    m_old       =    theta_old(1);                                 % number of segments
    h01_old     =    theta_old(2);                                 % cease-to-flow parameter of the first segment
    h_s_old     =    theta_old(2:m_old+2);                         % list fo break points (h01 is included in the list of break-points)
    h0_list_old =    [h01_old,theta_old(m_old+3:2*m_old+1)];       % list of cease-to-flow parameters
    a1_old      =    theta_old(2*m_old+2);                         % list of multiplier parameters
    b_list_old  =    theta_old(2*m_old+3:3*m_old+2);               % list of exponent parameters
    sigma2_old  =    theta_old(3*m_old+3:4*m_old+2);
    
    %% proposed parameters
    m_new       =    theta_new(1);
    h01_new     =    theta_new(2);
    h_s_new     =    theta_new(2:m_new+2);
    h0_list_new =    [h01_new,theta_new(m_new+3:2*m_new+1)];
    a1_new      =    theta_new(2*m_new+2);
    b_list_new  =    theta_new(2*m_new+3:3*m_new+2);
    sigma2_new  =    theta_new(3*m_new+3:4*m_new+2);
    
    %% density of m_new
    % probabilities (p(m_old)=0.8, p(m_old+1)=0.1,p(m_old-1)0.1)
    %
    m_max=floor(lambda*(hmax-hmin));
    m_min=1;
    if m_max~=m_min
        if m_new==m_old
            log_densm=log(0.8);
        elseif m_old==m_max || m_old==m_min
            log_densm=log(0.2);
        else
            log_densm=log(0.10);
        end
    else
        log_densm=0;
    end
    %}
    % if m is fixed
%     log_densm=0;
    %% density of h01_new
    if m_new==m_old
        % normal distribtuion
        %         half_range=min(abs(h01_old-h0_max),abs(h01_old-h0_min));
        %         sigma_h01=half_range/3;
        %         dens_h01=normpdf(h01_new-h01_old,sigma_h01);
        
        % beta distribution
        h01_min_tmp=max(h01_old-0.3,h0_min);
        h01_max_tmp=min(h01_old+0.3,h0_max);
        log_dens_h01 = betapdfAG((h01_new-h01_min_tmp)/(h01_max_tmp-h01_min_tmp),2,2);
    else
        log_dens_h01 = log(prior_h01_dens(h01_new,h0_min,h0_max));
    end
    
    %% density of h_s_new
    if m_new==m_old
        log_dens_h_s=0;
        for hs_ind=2:length(h_s_new)-1
            %             hs_tmp=h_s_old(hs_ind);
            %             hs_min=h_s_new(hs_ind-1);
            %             hs_max=h_s_old(hs_ind+1);
            %             half_range=min(abs(hs_tmp-hs_min),abs(hs_tmp-hs_max));
            %             sigma_hs=half_range/4;
            %
            %             dens_h_s=dens_h_s*normpdf(h_s_new(hs_ind)-hs_tmp,sigma_hs);
            
            % beta distribution
            hs_tmp=h_s_old(hs_ind);
            hs_min=max(hs_tmp-0.3,h_s_new(hs_ind-1));
            hs_max=min(hs_tmp+0.3,h_s_old(hs_ind+1));
            log_dens_h_s = log_dens_h_s + betapdfAG((h_s_new(hs_ind)-hs_min)/(hs_max-hs_min),2,2);
        end
    else
        log_dens_h_s = log(prior_h_s_dens(h_s_new,m_new,hmin,hmax));
    end
    
    %% density of h0j_new
    if m_new==m_old
        log_dens_h0j = 0;
        for h0j_ind = 2:m_new
            h0tmp = h0_list_old(h0j_ind);
            h0j_min = max(h0tmp-0.3,h0_min);
            h0j_max = min(h0tmp+0.3,h_s_new(h0j_ind));
            log_dens_h0j = log_dens_h0j + betapdfAG((h0_list_new(h0j_ind)-h0j_min)/(h0j_max-h0j_min),2,2);
        end
    else
        log_dens_h0j = log(prior_h0j_dens(h0_list_new(2:end),h0_min,h0_max,h_s_new));
    end
    %% density of exponent parameter
    if m_new==m_old
        log_densb=0;
        for b_ind=1:length(b_list_new)
            % normal distribution
            %             btmp=b_list_old(b_ind);
            %             half_range=min(abs(btmp-bmin),abs(btmp-bmax));
            %             sigma_b=half_range/3;
            %             densb=densb*normpdf(b_list_new(b_ind)-btmp,sigma_b);
            
            % beta distribution
            btmp=b_list_old(b_ind);
            bmin_tmp=max(btmp-0.3,bmin);
            bmax_tmp=min(btmp+0.3,bmax);
            log_densb = log_densb + betapdfAG((b_list_new(b_ind)-bmin_tmp)/(bmax_tmp-bmin_tmp),2,2);
        end
    else
        log_densb = log(prior_b_dens(b_list_new,bmin,bmax));
    end
    
    %% density of multiplier parameter
    if m_new==m_old
        % normal distribution
        %         half_range=min(abs(a1_old-amin),abs(a1_old-amax));
        %         sigma_a=half_range/3;
        %         densa1=normpdf(a1_new-a1_old,sigma_a);
        
        % beta distribution
        amin_tmp=max(a1_old-0.3,amin);
        amax_tmp=min(a1_old+0.3,amax);
        log_densa1 = betapdfAG((a1_new-amin_tmp)/(amax_tmp-amin_tmp),2,2);
        
    else
        log_densa1 = log(prior_a_dens(a1_new,amin,amax));
    end
    
    %% density of error variance parameter
    if m_new==m_old
        log_dens_sigma2 = 0;
        for sig_ind=1:m_new
            sigma2_max=sigma2_old(sig_ind)+10;
            sigma2_min=max(sigma2_old(sig_ind)-10,0);
            log_dens_sigma2 = log_dens_sigma2 + betapdfAG((sigma2_new(sig_ind)-sigma2_min)/(sigma2_max-sigma2_min),2,2);
        end
    else
        log_dens_sigma2 = log(prior_sigma2_dens(sigma2_new,alpha,beta));
    end
    %% log of proposed parameter set
    log_dens_prop_prior = log_densm + log_dens_h01 + log_dens_h_s + log_dens_h0j + log_densb + log_densa1 + log_dens_sigma2;
end