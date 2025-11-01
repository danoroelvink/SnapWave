close all
par={'hs','tp','wd','ds','zs'}
for i=1:length(par)
    a=load(['ERA5_',par{i},'.txt']);
    figure;
    plot(a(:,1)/3600/24,a(:,2:end))
    xlim([29 35])
    title(par{i})
end