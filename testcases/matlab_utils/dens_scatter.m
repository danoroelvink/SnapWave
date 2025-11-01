function [dens]=dens_scatter(valm,valc,n)
% DENS_SCATTER - assign a density to measured/computed data pairs

valmax=max(max(valm),max(valc));
dv=valmax/n;
count=zeros(n+1);
dens=zeros(size(valm));
for k=1:length(valm);
    i=max(round(valm(k)/dv),1);
    j=max(round(valc(k)/dv),1);
    count(i,j)=count(i,j)+1;
end
for k=1:length(valm);
    i=max(round(valm(k)/dv),1);
    j=max(round(valc(k)/dv),1);
    dens(k)=count(i,j);
end

end
