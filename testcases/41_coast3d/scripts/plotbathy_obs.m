clear all;close all
figure('units','normalized','position',[0.2 0.2 0.6 0.6],'color','w');
g(1).netname='..\run\DCSM\dcsm-fine-nl\UGdcsm-fm_100m2coast3d_v3_era5bnd_net.nc';
g(2).netname='..\run\LocalModel\UGegm51_net.nc';
hisname='..\output\dcsm-fine-nl\snapwave_his.nc';
xyb=load('..\run\DCSM\dcsm-fine-nl\ERA5_COAST3D_points.txt');
nstat=9;
for i=1:2
%% Read full net
g(i).xn=ncread(g(i).netname,'mesh2d_node_x');
g(i).yn=ncread(g(i).netname,'mesh2d_node_y');
g(i).zn=ncread(g(i).netname,'mesh2d_node_z');
if i==2
[g(i).xn,g(i).yn]=convertCoordinates(g(i).xn,g(i).yn,'CS2.code',4326,'CS1.code',28992);
end
%F=scatteredInterpolant(
face_nodes=ncread(g(i).netname,'mesh2d_face_nodes');
for k=1:size(face_nodes,2);
    if isnan(face_nodes(4,k));
        face_nodes(4,k)=face_nodes(1,k);
    end;
end
g(i).xface=g(i).xn(face_nodes);
g(i).yface=g(i).yn(face_nodes);
g(i).zface=g(i).zn(face_nodes);
subplot(1,2,i)
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
 p=patch(g(i).xface,g(i).yface,g(i).zface); hold on;
 %shading interp;
%scatter(xn,yn,5,zn,'filled');
%axis equal
xstat=ncread(hisname,'station_x')
ystat=ncread(hisname,'station_y')
namstat=ncread(hisname,'station_name')
xlabel('longitude (^o)')
ylabel('latitude (^o)')

xlim([4.602 4.622])
ylim([52.58,52.63])
caxis([-25 5])
c=colorbar
ylabel(c,'Bed level (m)','fontsize',9)

hold on
scatter(xstat,ystat,10,'r','filled')
%text(xstat,ystat,deblank(namstat'),'fontweight','bold','color','k',)
%scatter(xyb(:,1),xyb(:,2),20,'g','filled')
hBase = plot_openstreetmap('Alpha', 0.6, 'Scale', 2); % Basemap.
end

export_fig('..\results\Coast3D-detail.png','-dpng','-r300')
subplot(121)
xlim([4.608 4.614])
ylim([52.5975 52.6125])
%shading interp
nm=[1,4,5,6,7,8]
scatter(xstat(nm),ystat(nm),10,'r','filled')
t=text(xstat(nm),ystat(nm),deblank(namstat(:,nm)'),'fontsize',9,'fontweight','bold','color','k','rotation',45)
caxis([-6 0])
c=colorbar
ylabel(c,'Bed level (m)','fontsize',9)

subplot(122)
cla
export_fig('Coast3D-meas.png','-dpng','-r600')

%% Overall picture
if 1
figure('units','normalized','position',[0.2 0.2 0.6 0.6],'color','w');
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
p=patch(g(1).xface,g(1).yface,g(1).zface); hold on;
shading flat;
%axis equal
xlabel('longitude (^o)')
ylabel('latitude (^o)')

p=patch(g(2).xface,g(2).yface,g(2).zface); 
caxis([-25 5])
c=colorbar
ylabel(c,'Bed level (m)','fontsize',9)
%plot(xbox,ybox,'k','linewidth',1)
%scatter(xstat,ystat,5,'r','filled')
scatter(xyb(:,1),xyb(:,2),20,'g','filled')
xlim([min(g(1).xn),max(g(1).xn)])
ylim([min(g(1).yn),max(g(1).yn)])
hBase = plot_openstreetmap('Alpha', 0.4, 'Scale', 2); % Basemap.
export_fig('..\results\Coast3D-overview.png','-dpng','-r300')

end