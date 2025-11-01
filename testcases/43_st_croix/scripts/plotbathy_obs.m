clear all;close all
addpath('..\..\matlab_utils')
netname='..\run\UGera2cdip_stcroix_dryout_net.nc';
hisname='..\output\snapwave_dryout_his.nc';
xyb=load('..\run\ERA5_CDIP_points.txt');
nstat=2;
%% Read full net
xn=ncread(netname,'mesh2d_node_x');
yn=ncread(netname,'mesh2d_node_y');
zn=ncread(netname,'mesh2d_node_z');
[xn,yn]=convertCoordinates(xn,yn,'CS2.code',4326,'CS1.name','WGS 84 / UTM zone 20N');
%F=scatteredInterpolant(
face_nodes=ncread(netname,'mesh2d_face_nodes');
for i=1:size(face_nodes,2);
    if isnan(face_nodes(4,i));
        face_nodes(4,i)=face_nodes(1,i);
    end;
end
xface=xn(face_nodes);
yface=yn(face_nodes);
zface=zn(face_nodes);
zn(zn>0)=nan;
figure('units','normalized','position',[0.2 0.2 0.6 0.6],'color','w');
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
%p=patch(xface,yface,zface); 
shading interp;
scatter(xn,yn,10,zn);
hold on;
axis equal
% xstat=ncread(hisname,'station_x',1,6)
% ystat=ncread(hisname,'station_y',1,6)
xstat=ncread(hisname,'station_x');
ystat=ncread(hisname,'station_y');
%namstat=ncread(hisname,'station_name',[1,1],[6,6]);
namstat=[' Christiansted';' Fareham      ']
xlabel('longitude (^o)')
ylabel('latitude (^o)')
% xlabel('Easting (km)')
% ylabel('Northing (km)')
% 
caxis([-1000 100])
c=colorbar
ylabel(c,'Bed level (m)','fontsize',9)
[xstat,ystat]=convertCoordinates(xstat,ystat,'CS2.code',4326,'CS1.name','WGS 84 / UTM zone 20N');
hold on
scatter(xstat,ystat,40,'r','filled')
text(xstat,ystat,deblank(namstat),'fontsize',9,'color','r')
scatter(xyb(:,1),xyb(:,2),40,'g','filled')
hBase = plot_openstreetmap('Alpha', 0.4, 'Scale', 2); % Basemap.
%EHY_plot_satellite_map
xlim([-65 -64])
ylim([17.5 18])

print('..\results\StCroix_model_map.png','-dpng')

%% Overall picture
if 0
figure('units','normalized','position',[0.2 0.2 0.6 0.6],'color','w');
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
p=patch(xface,yface,zface); hold on;
shading flat;
axis equal
xlabel('longitude (^o)')
ylabel('latitude (^o)')

xbox=[5.47,5.75,5.75,5.47,5.47];
ybox=[53.43,53.43,53.55,53.55,53.43];
caxis([-25 5])
c=colorbar
ylabel(c,'Bed level (m)','fontsize',9)
plot(xbox,ybox,'k','linewidth',1)
scatter(xstat,ystat,5,'r','filled')
scatter(xyb(:,1),xyb(:,2),20,'g','filled')
xlim([min(xn),max(xn)])
ylim([min(yn),max(yn)])
hBase = plot_openstreetmap('Alpha', 0.4, 'Scale', 2); % Basemap.
end