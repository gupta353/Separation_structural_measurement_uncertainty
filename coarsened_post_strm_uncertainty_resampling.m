% this routine randomly selects 100000 samples from the distribution of
% streamflow samples obtained by the multivariate Gaussian distribution
% fitted to ML obtained streamflow series

clear all
close all
clc

direc = 'E:/covmat_Jefferys_prior';
alpha = 1;

% select the samples indices to be picked
ind = randsample(1:10^6,10^5);
ind = sort(ind);

% pick samples from first   
folder = ['strm_samples_alpha=',num2str(alpha)];
strm_samps = [];
for count =1:1000
    
    ind_range = (count-1)*1000+1:count*1000;
    cind = intersect(ind,ind_range);
    cind = cind - (count-1)*1000;
    
    % read the data containing the indices
    fname = [folder,'_',num2str(count),'.mat'];
    filename = fullfile(direc,folder,fname);
    load(filename);
    
    strm_samps = [strm_samps;strm_tmp(cind,:)];
    clear strm_tmp
end

sname = ['strm_samples_alpha=',num2str(alpha),'.mat'];
filename = fullfile(direc,sname);
save(filename)