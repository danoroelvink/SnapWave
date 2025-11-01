function [ bias, relbias, sci, rmse ] = skill( measured, computed )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% compute skills 

zmt     = measured; 
zct     = computed; 

% r2      = mean((zct-mean(zct)).*(zmt-mean(zmt)))/(std(zmt)*std(zct)); 
rmse    = sqrt(nanmean((zct-zmt).^2)); 
rmsm    = sqrt(nanmean(zmt.^2)); 
sci     = rmse/max(rmsm,abs(nanmean(zmt))); 

relbias = nanmean(zct-zmt)/max(rmsm,abs(nanmean(zmt))); 
bias    = nanmean(zct-zmt); 

end


  
 