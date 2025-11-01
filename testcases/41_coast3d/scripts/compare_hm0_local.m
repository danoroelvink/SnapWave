clear all;close all
addpath('..\..\matlab_utils');
load ('..\data\Coast3D_validationdata.mat')
fname='..\output\snapwave_his.nc';
s2=sqrt(2);
Hindex=[1,0,2,7,1,1,1,1,0];
Hfac=[1,0,s2,s2,s2,s2,s2,s2,0];
names=deblank(ncread(fname,'station_name')');
tc=ncread(fname,'time');
refstr=ncreadatt(fname,'time','units');
refstr=refstr(15:end)
reftime=datenum(refstr)
tc=reftime+double(tc/3600/24);
js=load('..\run\LocalModel\jonswap.txt');
tjs=reftime+js(:,1)/3600/24;
tstart=datenum(1998,11,1,0,0,0);
tjs(tjs<tstart)=nan;
tlast=tjs(end);
hm0_js=js(:,2);
%tlast=data(1).timeseries(Hindex(1)).time(end);
tc(tc>tlast)=nan;
hm0=ncread(fname,'point_hm0');
figure(2)
%% plot boundary conditions
% subplot(7,1,1);
% plot(tjs,hm0_js,'linewidth',2);
% hold on

i=1
s(1)=struct;
for ind=2:8
    if Hindex(ind)~=0&~strcmp(names(ind,:),'7a ')
        i=i+1;
        figure(2)
        set(gcf,'units','normalized','position',[0.2 0.2 0.6 0.6],'color','w');

        subplot(7,1,i)
        tm=data(ind).timeseries(Hindex(ind)).time;
        tm(tm>tlast)=nan;
        valm=data(ind).timeseries(Hindex(ind)).val*Hfac(ind);
        valm=valm(~isnan(tm));
        tm=tm(~isnan(tm));
        tm=tm(~isnan(valm));
        valm=valm(~isnan(valm));
        valc=squeeze(hm0(ind,:));
        valc(valc==-999)=nan;
        if ind>0
            plot(tc,valc,'b',tm,valm,'.r','linewidth',2);
        else
            plot(tm,valm,'.r');
        end
        title(names(ind,:));
        datetick
        grid on
        ylim([0 3])
        ylabel('Hm0 (m)')

        figure(3)
        set(gcf,'units','normalized','position',[0.2 0.2 0.6 0.6],'color','w');
        subplot(2,3,i)
        valcm=interp1(tc(~isnan(tc)),valc(~isnan(tc)),tm);
        bias(ind)=nanmean(valcm-valm(~isnan(valm)));
        rms(ind)=sqrt(nanmean((valcm-valm(~isnan(valm))).^2));
        meanm(ind)=nanmean(valm(~isnan(valm)));
        meanc(ind)=nanmean(valcm);
        hsm=valm;
        hsc=valcm;
        % hsm=hsm(~isnan(hsc));
        % hsc=hsc(~isnan(hsc));
        dens=dens_scatter(hsm,hsc,10);
        scatter(hsm,hsc,10,dens/max(dens),'filled');
        hold on
        plot([0 3],[0 3],'r','linewidth',1)
        set(gca,'color','k')
        [ bias, relbias, sci ] = skill(hsm,hsc);
        [S(i-1)]=skill0([tm,hsm],[tm,hsc]);
        pointname{i-1}=names(ind,:);
        text(0.5,2.5,['bias=' num2str(bias,'%.2f')],'FontSize',7,'color','w')
        text(0.5,2.2,['SCI=' num2str(sci,'%.2f')],'FontSize',7,'color','w')
        title(names(ind,:));
        if i>3
            xlabel('Hm0 obs. (m)')
        end
        if mod(i,3)==1
            ylabel('Hm0 comp. (m)')
        end
        
    end
end
figure(2)
export_fig('..\results\Coast3D_Hm0_timeseries_comparison.png')
figure(3)

%% Add colour legend
subplot(231);
c=colorbar;
ylabel(c,'relative point density')
set(gca,'visible','off')

export_fig('..\results\Coast3D_Hm0_scatterplot.png')

fi=fopen('..\results\Coast3D_stats_Hm0.txt','w');
fprintf(fi,'%6s %8s %8s %8s %8s %8s %8s\n',...
    'point','r2','sci','relbias','bias','skill','rmse')
    for ipt=1:length(S)
        fprintf(fi,'%6s %8.4f %8.3f %8.3f %8.3f %8.3f %8.3f\n', ...
            pointname{ipt}, ...
            S(ipt).r2,S(ipt).sci,S(ipt).relbias, ...
            S(ipt).bias,S(ipt).skill,S(ipt).rmse)
    end
fclose(fi)
