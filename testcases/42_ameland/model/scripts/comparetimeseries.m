clear all; close all
starttime=datenum(2008,11,5);
stoptime =datenum(2008,12,31);
datadir='..\..\data_ameland\';
resultdir='..\results\';
hisname='..\output\ameland_his.nc';
csv=dir([datadir,'*.csv']);
nstat=6;
for i=1:nstat;%length(csv)
    fname=[datadir,csv(i).name];
    m=importdata(fname,',',1);
    for it=1:size(m.data,1)
        m.time(it)=datenum(m.textdata(it+1,1));
    end
    intime=m.time>=starttime&m.time<stoptime;
    m.Hm0=m.data(:,1);
    m.Tm10=m.data(:,4);
    m.time=m.time(intime);
    m.Hm0=m.Hm0(intime);
    m.Tm10=m.Tm10(intime);
    m.Hm0(m.Hm0==-999)=nan;
    m.Tm10(m.Tm10==-999)=nan;
    data(i)=m;
end
figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 1],'color','w')
for i=1:length(data)
    subplot(nstat,1,i)
    plot(data(i).time,data(i).Hm0,'k.')
    hold on
    datetick('x','dd/mm')
    ylabel('Hm0 (m)')
    xtick=[733713,733720,733727,733734,733743,733750,733757,733764,733774];
    xticklabel=['01/11';'08/11';'15/11';'22/11';'01/12';'08/12';'15/12';'22/12';'01/01'];
    set(gca,'xtick',xtick,'xticklabel',xticklabel);
end
figure(2)
set(gcf,'units','normalized','outerposition',[0 0 1 1],'color','w')
for i=1:length(data)
    subplot(nstat,1,i)
    plot(data(i).time,data(i).Tm10*1.1,'k.')
    hold on
    datetick('x','dd/mm')
    xtick=[733713,733720,733727,733734,733743,733750,733757,733764,733774];
    xticklabel=['01/11';'08/11';'15/11';'22/11';'01/12';'08/12';'15/12';'22/12';'01/01'];
    set(gca,'xtick',xtick,'xticklabel',xticklabel)
    ylabel('Tm10 (s)')
end
fname=hisname;
time=nc_varget(fname,'time');
hm0=nc_varget(fname,'point_hm0');
tp=nc_varget(fname,'point_tp');
wavdir=nc_varget(fname,'point_wavdir');
names=nc_varget(fname,'station_name');
time=datenum(2008,9,2)+time/3600/24;
intime=time>=starttime&time<stoptime;
hm0=hm0(intime,:);
tp=tp(intime,:);
wavdir=wavdir(intime,:);
time=time(intime);
figure(1)
for i=1:nstat
    subplot(nstat,1,i)
    plot(time,hm0(:,i),'b','linewidth',2)
    title(deblank(names(i,:)))
end
figure(2)
for i=1:nstat
    subplot(nstat,1,i)
    plot(time,tp(:,i),'b','linewidth',2)
    title(deblank(names(i,:)))
end
%% Scatter plots and statistics
figure(3)
set(gcf,'units','normalized','outerposition',[0 0 1 1],'color','w')
t=tiledlayout(2,3)
for i=1:length(data)
    tm=data(i).time';
    Hm=data(i).Hm0;
    tc=time;
    Hc=hm0(:,i);
    dir=wavdir(:,1);
    Hcm=interp1(time,Hc,tm);
    dircm=interp1(time,dir,tm);
    tloc(i)=nexttile
    [xmap,ymap,zmap] = data_density_EQ(Hm,Hcm,50,50);
    [X,Y]=ndgrid(xmap,ymap);
    F=griddedInterpolant(X,Y,zmap);
    %scatter(Hm,Hcm,5,F(Hm,Hcm),'filled')
    scatter(Hm,Hcm,5,dircm,'filled')
    axis([0 max(Hm)*1.1 0 max(Hm)*1.1])
    set(gca,'color','k')
    hold on
    colormap parula;
    plot([0,1.1*max(Hm)],[0,1.1*max(Hm)],'w-','linewidth',1)
    xlabel('Hm0 obs. (m)')
    ylabel('Hm0 mod. (m)')
    RMS   = rms(Hcm-Hm);
    SI    = rms(Hcm-Hm)/rms(Hm);
    relbias  = nanmean(Hcm-Hm)/nanmean(Hm);
    skill = 1-nanvar(Hcm-Hm)/nanvar(Hm);
    title(deblank(names(i,:)))
    text(.05*max(Hm),.95*max(Hm),['SI        = ',num2str(SI,'%.2f')],'color','w','fontsize',9)
    text(.05*max(Hm),.88*max(Hm),['rel. bias = ',num2str(relbias,'%.2f')],'color','w','fontsize',9)
    text(.05*max(Hm),.81*max(Hm),['skill     = ',num2str(skill,'%.2f')],'color','w','fontsize',9)
end
c=colorbar(tloc(end))
ylabel(c,'wave direction (^oN)','FontSize',12)
figure(1);
export_fig([resultdir,'Ameland_timeseries_Hm0.png'],'-dpng')
figure(2);
export_fig([resultdir,'Ameland_timeseries_Tp.png'],'-dpng')
figure(3);
export_fig([resultdir,'Ameland_scatter_Hm0.png'],'-dpng')

for i=1:nstat
    S(i)=skill0([data(i).time',data(i).Hm0],[time,hm0(:,i)]);
end

fi=fopen([resultdir,'Ameland_stats_Hm0.txt'],'w');
fprintf(fi,'%6s %8s %8s %8s %8s %8s %8s\n',...
    'point','r2','sci','relbias','bias','skill','rmse')
for ipt=1:length(S)
    fprintf(fi,'%6s %8.4f %8.3f %8.3f %8.3f %8.3f %8.3f\n', ...
        deblank(names(ipt,:)), ...
        S(ipt).r2,S(ipt).sci,S(ipt).relbias, ...
        S(ipt).bias,S(ipt).skill,S(ipt).rmse)
end
fclose(fi)

for i=1:nstat
    S(i)=skill0([data(i).time',data(i).Tm10*1.1],[time,tp(:,i)]);
end

fi=fopen([resultdir,'Ameland_stats_Tp.txt'],'w');
fprintf(fi,'%6s %8s %8s %8s %8s %8s %8s\n',...
    'point','r2','sci','relbias','bias','skill','rmse')
for ipt=1:length(S)
    fprintf(fi,'%6s %8.4f %8.3f %8.3f %8.3f %8.3f %8.3f\n', ...
        deblank(names(ipt,:)), ...
        S(ipt).r2,S(ipt).sci,S(ipt).relbias, ...
        S(ipt).bias,S(ipt).skill,S(ipt).rmse)
end
fclose(fi)
