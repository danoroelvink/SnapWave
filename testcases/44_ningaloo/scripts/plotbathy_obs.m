clear all;close all
figure('units','normalized','position',[0.2 0.2 0.6 0.6],'color','w');
g(1).netname='..\run\ERA5_model\UGningaloo_5xref_net.nc';
g(2).netname='..\run\Local_model\UGningaloo_local_net.nc';
hisname='..\output\ningaloo_localfw60_his.nc';
xyb=load('..\run\ERA5_model\ERA5_bndpoints.txt');
nstat=9;
for i=1:2
%% Read full net
g(i).xn=ncread(g(i).netname,'mesh2d_node_x');
g(i).yn=ncread(g(i).netname,'mesh2d_node_y');
g(i).zn=ncread(g(i).netname,'mesh2d_node_z');
if 0
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

ylim([-22.31 -22.11])
xlim([113.8 113.9])
caxis([-25 5])
c=colorbar
ylabel(c,'Bed level (m)','fontsize',9)

hold on
scatter(xstat,ystat,10,'r','filled')
%text(xstat,ystat,deblank(namstat'),'fontsize',9,'color','r')
%scatter(xyb(:,1),xyb(:,2),20,'g','filled')
hBase = plot_openstreetmap('Alpha', 0.6, 'Scale', 3); % Basemap.
end

export_fig('..\results\Ningaloo-detail.png','-dpng','-r300')
subplot(131)
xlim([4.605 4.615])
ylim([52.602,52.606])
shading interp
scatter(xstat,ystat,10,'r','filled')
t=text(xstat,ystat,deblank(namstat'),'fontsize',9,'color','r','rotation',45)
caxis([-10 0])
c=colorbar
ylabel(c,'Bed level (m)','fontsize',9)

subplot(121)
cla

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
export_fig('..\results\Ningaloo-overview.png','-dpng','-r600')

end