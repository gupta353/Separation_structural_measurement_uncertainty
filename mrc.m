% This routine estimates the master recession curve
% 
% Ref: Lamb and Beven (1997), Tallaksen (1995)

clear all
close all
clc

direc='D:/Research/Thesis_work/Structural_vs_measurement_uncertainty/matlab_codes';
darea=2764.0831;                % drainage area in km2
rec_length_thresh=4;           % minimum lenght of recession period to be considered (in days)
prcp_strm_ratio_thresh=0.2;     % threshold for precipitation to streamflow ratio (recommended <=0.1)
evap_strm_ratio_thresh_pos=0.1;     % threshold for positive evaporation to streamflow ratio (recommended <=0.1)
evap_strm_ratio_thresh_neg=-0.1;     % threshold for negative evaporation to streamflow ratio (recommended <=0.1)
%% load streamflow, rainfall and evaporation data
% streamflow
strm_fname='streamflow_18.txt';
filename=fullfile(direc,'huc_04100003','streamflow_data',strm_fname);
fid=fopen(filename,'r');
data=textscan(fid,'%s%f','delimiter',char(32),'headerlines',1);
fclose(fid);
strm_date=data{1};
strm_vals=data{2}*0.028316847;      % cfs to cms
strm_vals=strm_vals/darea/1000*3600*24;     % cms to mm day^{-1}

% precipitation
prcp_fname='rainfall.txt';
filename=fullfile(direc,'huc_04100003',prcp_fname);
fid=fopen(filename,'r');
data=textscan(fid,'%s%s%f','delimiter',char(32),'headerlines',1);
fclose(fid);
prcp_date=data{1};
prcp_time=data{2};
prcp_vals=data{3};
% replace NaNs by means of previous and next day
indnan=find(isnan(prcp_vals));
for nan_ind=1:length(indnan)
    indtmp=[indnan(nan_ind)-1;indnan(nan_ind)+1];
    prcp_vals(indnan(nan_ind))=mean(prcp_vals(indtmp));
end
% convert rainfall into daily timescale from hourly timescale
prcp_tmp=reshape(prcp_vals,24,length(prcp_vals)/24);
prcp_vals=sum(prcp_tmp)';

% evaporation
evap_fname='evaporation.txt';
filename=fullfile(direc,'huc_04100003',evap_fname);
fid=fopen(filename,'r');
data=textscan(fid,'%s%f','delimiter',char(32),'headerlines',1);
fclose(fid);
evap_date=data{1};
evap_vals=data{2};

%% extarct recession periods
strm_shifted=strm_vals(1:end-1);
differ=strm_vals(2:end)-strm_shifted;
ind=find(differ<0);
rec_times=ind+1;

%
rec_time_shifted=rec_times(1:end-1);
differ=rec_times(2:end)-rec_time_shifted;
ind=find(differ>1);
rec_lengths=[ind(1);ind(2:end)-ind(1:end-1)];
break_times=rec_times(ind+1);

figure; hold on;
time_count=1;
for rec_ind=1:length(rec_lengths)
    
    rec_length_tmp=rec_lengths(rec_ind);
    plot_times=rec_times(time_count:time_count+rec_length_tmp-1);
    strm_vals_tmp=strm_vals(plot_times);
    
    rec_period{rec_ind}(:,1)=plot_times;
    rec_period{rec_ind}(:,2)=prcp_vals(plot_times);
    rec_period{rec_ind}(:,3)=evap_vals(plot_times);
    rec_period{rec_ind}(:,4)=strm_vals(plot_times);
    
    plot(plot_times,strm_vals_tmp);
    
%     pause(1);
    time_count=time_count+rec_length_tmp;
    
end
hold off;
xlabel('Time-steps (days)','fontname','arial','fontsize',12);
ylabel('Streamflow (mm day^{-1})','fontname','arial','fontsize',12);
title('Recession periods','fontname','arial','fontsize',12);
set(gca,'fontname','arial','fontsize',12);
%% Compute statistics to filter the recession periods
% compute total volumes of streamflow, rainfall, and evaporation during
% each recession period
for rec_ind=1:length(rec_period)
    
    prcp_vol(rec_ind)=sum(rec_period{rec_ind}(:,2));
    evap_vol(rec_ind)=sum(rec_period{rec_ind}(:,3));
    strm_vol(rec_ind)=sum(rec_period{rec_ind}(:,4));
        
end

prcp_strm_ratios=prcp_vol./strm_vol;
evap_strm_ratios=evap_vol./strm_vol;

% plot histograms of ratios
figure; hist(prcp_strm_ratios);
xlabel('Precipitation streamflow ratio','fontname','arial','fontsize',12);
ylabel('Number of samples in the bin','fontname','arial','fontsize',12);
set(gca,'fontname','arial','fontsize',12);

figure; hist(evap_strm_ratios);
xlabel('Evaporation streamflow ratio','fontname','arial','fontsize',12);
ylabel('Number of samples in the bin','fontname','arial','fontsize',12);
set(gca,'fontname','arial','fontsize',12);
%% filter the recession periods
ind_prcp=find(prcp_strm_ratios<=prcp_strm_ratio_thresh);        % precipitation threshold
ind_evap=find(evap_strm_ratios<=evap_strm_ratio_thresh_pos...
    & evap_strm_ratios>=evap_strm_ratio_thresh_neg);        % evaporation threshold
ind_length=find(rec_lengths>=rec_length_thresh);                % leength threshold
ind_final1=intersect(ind_prcp,ind_evap);
ind_final=intersect(ind_final1,ind_length);

filt_rec_periods=rec_period(ind_final);                         % filtered recession periods

%% create mrc based on filtered recession periods
% sort the filtered recession periods
figure; hold on
for rec_ind=1:length(filt_rec_periods)
    
    comp_val(rec_ind,1)=min(filt_rec_periods{rec_ind}(:,4));
    plot(filt_rec_periods{rec_ind}(:,4));
    
end
hold off;
xlabel('Arbitrary time-steps (days)','fontname','arial','fontsize',12);
ylabel('Streamflow (mm day^{-1})','fontname','arial','fontsize',12);
title('Filtered recession periods','fontname','arial','fontsize',12);
set(gca,'fontname','arial','fontsize',12);

comp_val=[comp_val,(1:length(comp_val))'];
comp_val=sortrows(comp_val);

% create MRC
mrcf=[];
for mrcf_ind=2:size(comp_val,1);
    
    % pick the lower and upper recession curve
    low_rec_ind=comp_val(mrcf_ind-1,2);
    up_rec_ind=comp_val(mrcf_ind,2);
    lower_rec=sort(filt_rec_periods{low_rec_ind}(:,4),'descend');
    upper_rec=sort(filt_rec_periods{up_rec_ind}(:,4),'descend');
    
    min_upper_rec=min(upper_rec);
    differ=abs(lower_rec-min_upper_rec);
    ind=find(differ==min(differ));
    mrcf=[lower_rec(ind:end);mrcf];
    
end
mrcf=[upper_rec;mrcf];
% mrcf=mrcf;
figure; plot(mrcf)
xlabel('Arbitrary time-steps (days)','fontname','arial','fontsize',12);
ylabel('Streamflow (mm day^{-1})','fontname','arial','fontsize',12);
title('Master recession curve','fontname','arial','fontsize',12);
set(gca,'fontname','arial','fontsize',12);

% save MRC
sname='MRC.txt';
filename=fullfile(direc,'huc_04100003','MRC',sname);
fid=fopen(filename,'w');
fprintf(fid,'%s\t%s\n','time-step(days)','Streamflow(cms)');
for write_ind=1:length(mrcf)
    
    fprintf(fid,'%f\t%f\n',write_ind-1,mrcf(write_ind)*darea*1000/24/3600); % coversion from mm day^{-1} to m^3 s^{-1}
    
end
close(fid)
