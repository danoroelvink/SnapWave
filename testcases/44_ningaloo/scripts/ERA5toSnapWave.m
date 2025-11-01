clear all;close all
load('..\XBeach_data\XBEACH_data_v2.mat');
starttime=datenum(2009,6,12,0,0,0);
endtime=datenum(2009,6,30,23,0,0);
time=time-8/24;
ok=time>=starttime&time<==endtime;
subplot(511)
plot(time,Hsig_swell,'k');datetick; 
hold on
ylabel('Hs (m)')
subplot(512)
plot(time,tp,'k');datetick
ylabel('Tp (s)')
subplot(513)
plot(time,wd);datetick
ylabel('dir (^oN)')
ylim([180 360])
subplot(514)
plot(time,ds);datetick
ylabel('dir. spr. (^o)')
out=[tsec,hs]; save('ERA5_hs.txt','out','-ascii');
out=[tsec,tp]; save('ERA5_tp.txt','out','-ascii');
out=[tsec,wd]; save('ERA5_wd.txt','out','-ascii');
out=[tsec,ds]; save('ERA5_ds.txt','out','-ascii');
%% GTSM water level data
fname='ERA5\reanalysis_waterlevel_hourly_2009_06_v1.nc';
xgtsm=nc_varget(fname,'station_x_coordinate');
ygtsm=nc_varget(fname,'station_y_coordinate');
for i=1:length(lon);
    [mindist,indgtsm(i)]=min(hypot(xgtsm-lon(i),ygtsm-lat(i)));
end
lonloc = 114;
latloc = -22;
[mindist,indloc]=min(hypot(xgtsm-lonloc,ygtsm-latloc));

month={'06'};
zsgtsm=[];
zsloc=[];
timegtsm=[];
for i=1:length(month)
    fname=['ERA5\reanalysis_waterlevel_hourly_2009_',month{i},'_v1.nc'];
    zs=nc_varget(fname,'waterlevel');
    zsl=zs(:,indloc)
    zs=zs(:,indgtsm);
    time=nc_varget(fname,'time');
    time=datenum(1900,1,1)+time/3600/24;
    zsgtsm=[zsgtsm;zs];
    zsloc=[zsloc;zsl];
    timegtsm=[timegtsm;time];
end
%% There may be NaNs in these timeseries; remove.
for p=1:size(zsgtsm,2)
    for t=2:size(zsgtsm,1)-1
        if isnan(zsgtsm(t,p))
            zsgtsm(t,p)=.5*(zsgtsm(t-1,p)+zsgtsm(t+1,p));
        end
    end
end

intime=timegtsm>=starttime&timegtsm<endtime;
timegtsm=timegtsm(intime);
zsgtsm=zsgtsm(intime,:);
zsloc=zsloc(intime);
tsec=(timegtsm-timegtsm(1))*3600*24;
out=[tsec,zsgtsm];
save('ERA5_zs.txt','out','-ascii')
subplot(515)
plot(timegtsm,zsgtsm,timegtsm,zsloc);datetick
legend('ERA','local')
ylabel('zs (m)')
%% copy local time series to all boundaries
for i=1:size(zsgtsm,2)
    zsgtsm(:,i)=zsloc;
end
out=[tsec,zsgtsm];
save('GTSM_Ningaloo_zs.txt','out','-ascii')

