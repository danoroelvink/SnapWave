clear all;close all
load('..\XBeachdata\XBEACH_data_v2.mat');
starttime=datenum(2009,6,12,0,0,0);
endtime=datenum(2009,6,30,23,0,0);
time=time-8/24;
ok=time>=starttime&time<=endtime&~isnan(Tp_swell);
time=time(ok)';
tsec=(time-starttime)*24*3600;
Hs=Hsig_swell(1,:);
Hs=Hs(ok)';
Tp=Tp_swell(ok)';
dir=DirTp_S(ok)';
dep=depth(1,ok)';
time=time-datenum(2009,1,1);
dur=endtime-starttime;
window=round(dur-mod((dur)*24,25));
meandepth=mean(dep(1:window));
level=dep-meandepth;


figure;
subplot(411)
plot(time,Hs,'k','linewidth',1);
xlim([160 180])
hold on
ylabel('Hs (m)')
subplot(412)
plot(time,Tp,'k','linewidth',1);
xlim([160 180])
ylabel('Tp (s)')
subplot(413)
plot(time,dir,'k','linewidth',1);
xlim([160 180])
ylabel('dir (^oN)')
ylim([220 290])
subplot(414)
plot(time,level,'k','linewidth',1);
xlim([160 180])
ylabel('water level (m)')
xlabel('yearday')
ds=zeros(size(Hs))+20;
out=[tsec,Hs,Hs]; save('ningaloo_hs.txt','out','-ascii');
out=[tsec,Tp,Tp]; save('ningaloo_tp.txt','out','-ascii');
out=[tsec,dir,dir]; save('ningaloo_wd.txt','out','-ascii');
out=[tsec,ds,ds]; save('ningaloo_ds.txt','out','-ascii');
out=[tsec,level,level]; save('ningaloo_zs.txt','out','-ascii');
