
clear all;close all
load('..\output\f32har01.mat');
snapname='..\output\hv_map.nc';
Hm0=ncread(snapname,'hm0');
Hm0=Hm0(:,1);
dHm0=Hm0-Hsig;
xn=ncread(snapname,'mesh2d_node_x');
yn=ncread(snapname,'mesh2d_node_y');
zn=ncread(snapname,'mesh2d_node_z');
zn(zn>0)=0;

F=scatteredInterpolant(double(Xp(:)),double(Yp(:)),double(Hsig(:)));
Hsig=F(xn,yn);
dHm0=Hm0-Hsig;

face_nodes=ncread(snapname,'mesh2d_face_nodes');
for k=1:size(face_nodes,2);
    if isnan(face_nodes(4,k));
        face_nodes(4,k)=face_nodes(1,k);
    end
end
xface=xn(face_nodes)/1000;
yface=yn(face_nodes)/1000;
zface=zn(face_nodes);
Hm0f=Hm0(face_nodes);
Hsigf=Hsig(face_nodes);
dHm0f=dHm0(face_nodes);


figure;

subplot(222)
patch(xface,yface,Hm0f);shading interp;
colormap parula
c=colorbar;
ylabel(c,'Hm0 (m)')
grid on
axis equal
title('Hm0 SnapWave')
xlim([7 23])
ylim([2 20])

subplot(221)
patch(xface,yface,Hsigf);shading interp;
colormap parula
c=colorbar;
ylabel(c,'Hm0 (m)')
ylabel('Northing (km)')
grid on
axis equal
title('Hm0 unSWAN')
xlim([7 23])
ylim([2 20])

subplot(224)
patch(xface,yface,dHm0f);shading interp;
% hold on
% xl=[0,17,17];yl=[18,18,0];
% plot(xl,yl,'k','linewidth',1)
colormap parula
caxis([-.2 .2])
c=colorbar;
ylabel(c,'\Delta Hm0 (m)')
xlabel('Easting (km)')
title('\Delta Hm0 SnapWave-unSWAN')
axis equal
xlim([7 23])
ylim([2 20])
subplot(223);
patch(xface,yface,zface);shading interp;
axis equal;
caxis([-10 0])
c=colorbar
ylabel(c,'Bed level (m)')
ylabel('Northing (km)')
xlabel('Easting (km)')
title('Bed level')
xlim([7 23])
ylim([2 20])

print('..\results\SnapWave_vs_unSWAN_planview.png','-dpng','-r600')

figure;
ok=zn>-20&yn<18000&yn>5000&Hm0>0&Hm0<3.2;
scatter(Hsig(ok),Hm0(ok),5,zn(ok),'filled');
colormap parula
hold on;
plot(Hsig,Hsig,'k','linewidth',1)
xlabel('Hm0 unSWAN (m)')
ylabel('Hm0 SnapWave (m)')
grid on
c=colorbar
ylabel(c,'Depth (m)')

Hcm=Hm0(ok);
Hm=(Hsig(ok));
RMS   = rms(Hcm-Hm);
SI    = rms(Hcm-Hm)/rms(Hm);
relbias  = nanmean(Hcm-Hm)/nanmean(Hm);
skill = 1-nanvar(Hcm-Hm)/nanvar(Hm);
text(.05*max(Hm),.95*max(Hm),['SI        = ',num2str(SI,'%.2f')],'color','k','fontsize',9)
text(.05*max(Hm),.88*max(Hm),['rel. bias = ',num2str(relbias,'%.2f')],'color','k','fontsize',9)
text(.05*max(Hm),.81*max(Hm),['skill     = ',num2str(skill,'%.2f')],'color','k','fontsize',9)
print('..\results\SnapWave_vs_unSWAN_scatter.png','-dpng','-r600')

