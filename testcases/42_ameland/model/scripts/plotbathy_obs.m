clear all;
netname='..\run\UGameland_bathy2008_casulli_net.nc';
hisname='..\output\ameland_his.nc';
xyb=load('..\run\ERA5_bndpoints.txt');
nstat=6;
%% Read full net
xn=ncread(netname,'mesh2d_node_x');
yn=ncread(netname,'mesh2d_node_y');
zn=ncread(netname,'mesh2d_node_z');
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
figure('units','normalized','position',[0.2 0.2 0.6 0.6],'color','w');
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
p=patch(xface,yface,zface); hold on;
shading flat;
axis equal
xstat=ncread(hisname,'station_x',1,6)
ystat=ncread(hisname,'station_y',1,6)
namstat=ncread(hisname,'station_name',[1,1],[6,6])
xlabel('longitude (^o)')
ylabel('latitude (^o)')

xlim([5.47,5.75])
ylim([53.43,53.55])
caxis([-25 5])
c=colorbar
ylabel(c,'Bed level (m)','fontsize',9)

hold on
scatter(xstat,ystat,20,'r','filled')
text(xstat,ystat,deblank(namstat'),'fontsize',9,'color','r')
hBase = plot_openstreetmap('Alpha', 0.4, 'Scale', 2); % Basemap.
print('..\results\ameland_detail.png','-dpng','-r600')

%% Overall picture
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
print('..\results\ameland_overall.png','-dpng','-r600')