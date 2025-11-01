%% Evaluate model performance StCroix model, Snapwave version
% version JR
clear all;  close all;fclose all;
scenario  = {'CDIP_check'};%,...
scenariolegend  = {'observed','Computed'};%,...
col={'k','b','g'}
datadir   = '..\data\CDIP';
plotme=0;

% 1. Time administration
refdatestr   = '20100601000000';
stopdatestr  = '20101101000000';
startdate    = datenum(refdatestr,'yyyymmddHHMMSS');
stopdate     = datenum(stopdatestr,'yyyymmddHHMMSS') ;

% 2. Get data in model period

datafiles = {'christiansted','fareham'};

for ist = 1:length(datafiles)
    raw=importdata([fullfile(datadir,datafiles{ist}),'.txt'],' ',2);
    raw=raw.data;
    data(ist).datetime=datenum([raw(:,1:5),0*raw(:,5)]);
    data(ist).hs=raw(:,6);
    data(ist).tp=raw(:,7);
    data(ist).wd=raw(:,8);
end


% 3a. Get coordinates
ncfil  = '..\output\snapwave_his.nc';
X = nc_varget(ncfil,'station_x');
x  = X(:);
Y = nc_varget(ncfil,'station_y');
y  = Y(:);

% 3b. find the closest model point to each observation station
nobs = length(x);
idp  = zeros(nobs,1);

for id = 1:nobs
    % 3c. Retrieve model output in points
    nc     = fullfile(ncfil);
    time   = nc_varget(nc,'time');
    timexb = startdate + (time)./24./3600 ;
    
    % hrms
    tmp = nc_varget(nc,'point_hm0');  % time points
    hm0(:,id) = tmp(:,id);
    
    % Tp
    tmp = nc_varget(nc,'point_tp');  % time points
    tp(:,id) = tmp(:,id);
    
    % Dir
    tmp = nc_varget(nc,'point_wavdir');  % time points
    wd(:,id) =  tmp(:,id);
end

%% 4b. Plot time series
% scrsz = get(groot,'ScreenSize');
figure;
set(gcf, 'DefaultLineLineWidth', 1,'Color','w');
set(gcf,'units','normalized','outerposition',[0 0.2 1 .6])
for irun=1:length(scenario)
    for id = 1:nobs
        
        a1 = subplot(2,1,id);
        sel=data(id).datetime>=startdate&data(id).datetime<=stopdate&~isnan(data(id).hs);
        datem=data(id).datetime(sel);
        hsm=data(id).hs(sel)
        plot(datem,hsm,'.b',timexb,hm0(:,id),col{irun});
        SK(id)=skill0([datem,hsm],[timexb,hm0(:,id)]);
        datetick
        ylabel('Hm0 (m)')
        %ylim([0 1.5])
        hold on
        if irun==length(scenario)&id==nobs
            legend(scenariolegend)
        end
    end
end
fname='..\results\StCroix_validation_timeseries.png';

export_fig(fname,'-dpng','-r500');
%% 5b. Plot scatter plots

figure
set(gcf, 'DefaultLineLineWidth', 1,'Color','w');
for irun=1:length(scenario)
    for id = 1:nobs
        
        a1 = subplot(2,nobs,2*id-1);hold on
        sel=data(id).datetime>=startdate&data(id).datetime<=stopdate&~isnan(data(id).hs);
        datem=data(id).datetime(sel);
        hsm=data(id).hs(sel)
        hsc=interp1(timexb,hm0(:,id),datem)
        %scatter(hsm,hsc,10,'filled','CData',col(i5,:),'DisplayName','observed');
        %plot(hsm,hsc,'.','color',col{irun});
        dens=dens_scatter(hsm,hsc,100);
        scatter(hsm,hsc,5,dens,'filled')
        plot([0 4.5],[0 4.5],'r','linewidth',1)
        set(gca,'color','k')
        [ bias, relbias, sci ] = skill(hsm,hsc);
        text(0.5,4,['bias=' num2str(bias,'%.2f')],'FontSize',7,'color','w')
        text(0.5,3.5,['SCI=' num2str(sci,'%.2f')],'FontSize',7,'color','w')
        if id == 1;ylabel('Hm0 [m]');end
        set(gca,'XTickLabel',[]);
        if id <= 8
            title(sprintf('Station %s', datafiles{id}))
        end
        %xlim([0 1]);ylim([0 1]);
        box on
        
        
        a3 = subplot(2,nobs,2*id);hold on
        tpm=data(id).tp(sel)
        tpc=interp1(timexb,tp(:,id),datem)
        %scatter(tpm,tpc,10,'filled','CData',col(i5,:),'DisplayName','observed');
        %plot(tpm,tpc,'.','color',col{irun});
        dens=dens_scatter(tpm,tpc,100);
        scatter(tpm,tpc,5,dens,'filled')
        set(gca,'color','k')
        plot([0 20],[0 20],'r','linewidth',1)
        [ bias, relbias, sci ] = skill( tpm,tpc );
        text(2,16,['bias=' num2str(bias,'%.2f')],'FontSize',7,'color','w')
        text(2,14,['SCI=' num2str(sci,'%.2f')],'FontSize',7,'color','w')
        
        if id == 1; ylabel('Tp [s]');end
        xlabel('observed [m]');
        xlim([0 20]); ylim([0 20]);box on
        if id == 8;legend;end
    end
end
fname='..\results\StCroix_scatterplot.png';

export_fig(fname,'-dpng','-r500');

fi=fopen(['..\results\StCroix_stats_Hm0.txt'],'w');
fprintf(fi,'%20s %8s %8s %8s %8s %8s %8s\n',...
    'point','r2','sci','relbias','bias','skill','rmse')
    for ipt=1:length(SK)
        fprintf(fi,'%20s %8.4f %8.3f %8.3f %8.3f %8.3f %8.3f\n', ...
            deblank(datafiles{ipt}), ...
            SK(ipt).r2,SK(ipt).sci,SK(ipt).relbias, ...
            SK(ipt).bias,SK(ipt).skill,SK(ipt).rmse)
    end
fclose(fi)

function [dens]=dens_scatter(valm,valc,n)
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
