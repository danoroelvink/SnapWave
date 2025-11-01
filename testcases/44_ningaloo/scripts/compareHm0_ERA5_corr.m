clear all;close all;
addpath('..\..\matlab_utils');
load('..\data\spectral_elevation_results_latestRyan.mat');
col={'b','r','.k'};
label={'A1','B1','C1','C1','C3','C4','C5','C6','E6'}
fname='..\output\ningaloo_fw60_his.nc';
%sp=[5,6,3,3,7,11,15,19,8];
sp=[0,0,1,1,2,3,4,5,0];
figure(1)
for i=1
    names=nc_varget(fname,'station_name');
    hm0=ncread(fname,'point_hm0');
    time=double(ncread(fname,'time'));
    units=ncreadatt(fname,'time','units');
    reftime=datenum(units(15:24));
    time=reftime+time/3600/24;
    day_model=time-datenum(2009,1,1)
    for j=3:8%1:size(d(i).hm0,1)
        subplot(5,1,sp(j))
        plot(day_model,hm0(j,:),col{i},'linewidth',2);
        hold on
    end

    %day_obs=time_burst-8/24-datenum(2009,1,1);
    day_obs=time_burst-7/24-datenum(2009,1,1);
end
hsall_obs=[];
hsall_mod=[];
i=0;
for j=3:8;%1:size(d(i).hm0,1)
    %subplot(5,4,sp(j))
    figure(1)
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'color','w')

    subplot(5,1,sp(j))
    k=j;
    if j==9
        k=j+1
    end
    plot(day_obs,Hsig_swell(k,:),col{3},'linewidth',2);
    grid on
    ylabel('Hm0 (m)')
    if j==8
        xlabel('yearday')
    end
    xlim([160 180])
    yl=get(gca,'ylim');
    text(161,.8*yl(2),label{j})
    hsc=interp1(day_model,hm0(j,:),day_obs)
    hsm=Hsig_swell(k,:);
    figure(2)
    set(gcf,'units','normalized','outerposition',[0 0 1 1],'color','w')

    if j~=4
    subplot(2,3,sp(j))
    dens=dens_scatter(hsm,hsc,20);
    scatter(hsm,hsc,10,dens,'filled')
    hold on
    yl=[0 max(max(hsc),max(hsm))]
    plot([0 3],[0 3],'r','linewidth',1)
    xlim(yl);
    ylim(yl);
    set(gca,'color','k')
    [ bias, relbias, sci ] = skill(hsm,hsc);
    i=i+1;
    SK(i)=skill0([time_burst',hsm'],[time_burst',hsc']);
    name{i}=label{j};
    text(0.1*yl(2), .8*yl(2),['relbias=' num2str(relbias,'%.2f')],'FontSize',7,'color','w')
    text(0.1*yl(2), .7*yl(2),['SCI=' num2str(sci,'%.2f')],'FontSize',7,'color','w')
    title(label{j});
    end

end

figure(1)
export_fig('..\results\ningaloo_timeseries_C_series_ERA5_corr.png','-dpng','-r300')
figure(2)
export_fig('..\results\ningaloo_scatter_C_series_ERA5_corr.png','-dpng','-r300')


fi=fopen(['..\results\Ningaloo_ERA5_stats_Hm0_corr.txt'],'w');
fprintf(fi,'%6s %8s %8s %8s %8s %8s %8s\n',...
    'point','r2','sci','relbias','bias','skill','rmse')
    for ipt=1:length(SK)
        fprintf(fi,'%6s %8.4f %8.3f %8.3f %8.3f %8.3f %8.3f\n', ...
            deblank(name{ipt}), ...
            SK(ipt).r2,SK(ipt).sci,SK(ipt).relbias, ...
            SK(ipt).bias,SK(ipt).skill,SK(ipt).rmse)
    end
fclose(fi)

g=load('..\run\ERA5_model\ERA5_zs.txt');
l=load('..\run\Local_model\ningaloo_zs.txt');
tg=g(:,1)/3600;
tl=l(:,1)/3600;
zsg=g(:,2);
zsl=l(:,2);
figure;
set(gcf,'color','w')
subplot(211)
plot(tg-1,zsg,tl,zsl)
title ('corrected wl')
xlim([0 450])
xlabel('time (hr)')
ylabel('zs (m)')
subplot(212)
plot(tg,zsg,tl,zsl)
title ('uncorrected wl')
xlim([0 450])
xlabel('time (hr)')
ylabel('zs (m)')
export_fig('..\results\zs_comparison_Ningaloo.png','-dpng')



