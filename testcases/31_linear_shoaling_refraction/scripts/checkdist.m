crit=20;
fname='shoalref_fine_neu_map.nc';
ee=ncread(fname,'ee');
x=ncread(fname,'mesh2d_node_x');
y=ncread(fname,'mesh2d_node_y');
eec=squeeze(ee(:,abs(y-5000)<crit,:));
eec3=squeeze(eec(:,:,3));
xc=x(abs(y-5000)<crit);
 theta=[1:180]*1-0.5-90+60;
figure;pcolor(xc,theta,eec3)
hold on;
plot(xhis,dirm(3,:),'k','linewidth',2)
figure;
plot(theta,eec3)
