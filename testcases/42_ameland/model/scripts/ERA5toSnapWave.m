clear all
fname='..\ERA5\ERA5_NL_2008.nc';
starttime=datenum(2008,9,1,0,0,0);
endtime=datenum(2008,12,31,23,0,0);
%% read zs from jonswap.txt file
% jons=load('..\jonswap.txt');
% tsecjons=jons(:,1);
% tjons=starttime+tsecjons/3600/24;
% zsjons=jons(:,end);
%fname='adaptor.mars.internal-1678378664.7782328-19880-7-f0cfd7d0-9495-445f-b5ca-afcbc8964a4e.nc';
%fname='adaptor.mars.internal-1678701598.9426239-1186-18-899c6f73-e554-4019-abe9-0668570e40d5.nc'
nc_dump(fname)
xyp=load("ERA5_bndpoints.txt");
lon=xyp(:,1);
lat=xyp(:,2);
x=nc_varget(fname,'longitude')
y=nc_varget(fname,'latitude')
for i=1:length(lon);indlon(i)=find(x==lon(i));end
for i=1:length(lat);indlat(i)=find(y==lat(i));end
swh=nc_varget(fname,'swh');
pp1d=nc_varget(fname,'pp1d');
mwd=nc_varget(fname,'mwd');
wdw=nc_varget(fname,'wdw');
for i=1:length(indlon)
    hs(:,i)=swh(:,indlat(i),indlon(i));
    tp(:,i)=pp1d(:,indlat(i),indlon(i));
    wd(:,i)=mwd(:,indlat(i),indlon(i));
    ds(:,i)=wdw(:,indlat(i),indlon(i))*180/pi;
end
time=nc_varget(fname,'time');
time=datenum(1900,1,1)+time/24;

intime=time>=starttime&time<endtime;
time=time(intime);
hs=hs(intime,:);
tp=tp(intime,:);
wd=wd(intime,:);
ds=ds(intime,:);
tsec=(time-time(1))*3600*24;
figure;
subplot(511)
plot(time,hs);datetick
ylabel('Hs (m)')
subplot(512)
plot(time,tp,'k');datetick
ylabel('Tp (s)')
subplot(513)
plot(time,wd);datetick
ylabel('dir (^oN)')
subplot(514)
plot(time,ds);datetick
ylabel('dir. spr. (^o)')
out=[tsec,hs]; save('ERA5_hs.txt','out','-ascii');
out=[tsec,tp]; save('ERA5_tp.txt','out','-ascii');
out=[tsec,wd]; save('ERA5_wd.txt','out','-ascii');
out=[tsec,ds]; save('ERA5_ds.txt','out','-ascii');
%% GTSM water level data
fname='..\ERA5\reanalysis_waterlevel_hourly_2008_09_v1.nc';
xgtsm=nc_varget(fname,'station_x_coordinate');
ygtsm=nc_varget(fname,'station_y_coordinate');
for i=1:length(lon);
    [mindist,indgtsm(i)]=min(hypot(xgtsm-lon(i),ygtsm-lat(i)));
end
lonloc = 5.477685;
latloc = 53.530063;
[mindist,indloc]=min(hypot(xgtsm-lonloc,ygtsm-latloc));

month={'09','10','11','12'};
zsgtsm=[];
zsloc=[];
timegtsm=[];
for i=1:length(month)
    fname=['..\ERA5\reanalysis_waterlevel_hourly_2008_',month{i},'_v1.nc'];
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

tsec=(timegtsm-timegtsm(1))*3600*24;
out=[tsec,zsgtsm];
save('ERA5_zs.txt','out','-ascii')
subplot(515)
plot(timegtsm,zsgtsm(:,5),timegtsm,zsloc);datetick
legend('ERA','local')
ylabel('zs (m)')
%% copy local time series Ameland to all boundaries
for i=1:size(zsgtsm,2)
    zsgtsm(:,i)=zsloc;
end
out=[tsec,zsgtsm];
save('GTSM_Ameland_zs.txt','out','-ascii')

