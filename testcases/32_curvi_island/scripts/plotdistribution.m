clear all;close all;
fname='..\output\snapwave_map.nc';
x0=[-200,0,200,0];y0=[0,200,0,-200];
ee=ncread(fname,'ee');
H=ncread(fname,'hm0');
theta=ncread(fname,'theta');
ntheta=size(theta,1);
x=ncread(fname,'mesh2d_node_x');
y=ncread(fname,'mesh2d_node_y');
z=ncread(fname,'mesh2d_node_z');
fn=ncread(fname,'mesh2d_face_nodes');
for i=1:size(fn,2);
    if isnan(fn(4,i));
        fn(4,i)=fn(1,i);
    end;
end
H=H(:,5);%size(H,2));
[Hmax,imax]=max(H);
xn=x(fn);
yn=y(fn);
zn=z(fn);
Hn=H(fn);
figure('units','normalized','position',[0 0 1 .5],'color','w');
%scatter(x,y,5,squeeze(ee(ntheta/2,:,6)));axis equal;colorbar
%scatter(x,y,10,H(:,6),'filled');axis equal;colorbar
% colormap jet
% caxis([-max(H(:,6)) max(H(:,6))])
subplot(121)
p=patch(xn,yn,zn); hold on;
%plot(x(imax),y(imax),'ok')
colorbar
axis equal
ylim([-500 500])
title('bed level (m)')
xlabel('x (m)')
ylabel('y (m)')
subplot(122)
p=patch(xn,yn,Hn); hold on;
plot(x(imax),y(imax),'ok')
shading interp;
colorbar
axis equal
ylim([-500 500])
xlabel('x (m)')
[mindist,imin]=min(hypot(x-x0,y-y0));
x2=x;
x2(imin)=1e10;
[mindist2,imin2]=min(hypot(x2-x0,y-y0));
hold on;
plot(.5*(x(imin)+x(imin2)),.5*(y(imin)+y(imin2)),'ro','linewidth',2)
text(x0'+15,y0',num2str([1:4]'))
title('Hm0 (m)')
print('..\results\curvi_planview.png','-dpng','-r300')

figure
for j=5:-1:1

for i=1:4
    subplot(2,2,i);
    [~,ind]=sort(theta(:,j));
    plot(theta(ind,j),squeeze(.5*(ee(ind,imin(i),j)+ee(ind,imin2(i),j))),'linewidth',2);
    title(['point ',num2str(i)])
    xlim([0 360]);
    set(gca,'xtick',[0 090 180 270 360])
    grid on
    hold on
    if i==2
        legend('sweep 5','sweep 4','sweep 3','sweep 2','sweep 1','location','northwest')
    end
    xlabel('direction (^oN)')
    ylabel('ee (J/m2/rad)')
end
end
print('..\results\curvi_sweeps.png','-dpng','-r300')
