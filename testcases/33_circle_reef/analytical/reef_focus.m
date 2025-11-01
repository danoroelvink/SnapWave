clear all;close all;
R=350;d0=100;d1=1.5;T=12;dy=5;dx=5;dtheta=1;N=100;critdist=4.5*dx;
thetabin=2;thetabin0=-86;nthetabin=-2*thetabin0/thetabin;
g=9.81;thr=0.05;
xpt=[-90,0,90,0];ypt=[0,90,0,-90];
Nx=2*R/dx+1+200/dx;
Ny=2*R/dy+1+200/dy;
omega=2*pi/T;
[k0,C0,Cg0]=disper(d0,T);
[k1,C1,Cg1]=disper(d1,T);
y0=[-R-100:dy:R+100];
out=abs(y0)>R;
%% angle between ray and reef
theta0=asind(y0/R);
theta0(out)=90;
%% x coordinate where ray crosses reef
x0=-R.*cosd(theta0);
%% refracted wave direction relative to reef edge normal
thetar=asind(sind(theta0)*k0/k1);
thetar(out)=theta0(out);
%% Shoaling and refraction factor energy
Efac=Cg0*cosd(theta0)./(Cg1*cosd(thetar));
Efac(out)=1;
%% refracted wave direction relative to x axis
theta1=thetar-theta0;
%% line projected from reef edge
x1=x0+2.1*R*cosd(theta1);
y1=y0+2.1*R*sind(theta1);
%% coordinates of reef edge
theta=[0:dtheta:360];
xR=R*cosd(theta);
yR=R*sind(theta);
%% cross x1,y1 with reef edge
for m=1:length(x0)
    [tempx,tempy,~,~,~,indj,~,uj]=get_intersections([x0(m) x1(m)],[y0(m) y1(m)],xR,yR);
    if ~isempty(tempx)
       [x2(m),imax]=max(tempx);
       y2(m)=tempy(imax);
       thetan2=(1-uj(imax))*theta(indj(imax))+uj(imax)*theta(indj(imax)+1);
       theta2=atan2d(y2(m)-y0(m),x2(m)-x0(m));
       dtheta2=theta2-thetan2;
       dtheta3=asind(sind(dtheta2)*k1/k0);
       theta3(m)=dtheta3+thetan2;
       x3(m)=real(x2(m)+2*R*cosd(theta3(m)));
       y3(m)=real(y2(m)+2*R*sind(theta3(m)));
    else
       x2(m)=R;
       y2(m)=y0(m);
       x3(m)=3*R;
       y3(m)=y0(m);
    end
end
%% Plotting
figure('units','normalized','position',[0 0 1 .7],'color','w');
subplot(121)
plot(x0,y0,'.k')
axis equal
hold on;
plot([0*x0-450;x0],[y0;y0],'r');
plot([x0;x2],[y0;y2],'k');
plot([x2;x3],[y2;y3],'b');
plot(xR,yR,'linewidth',2)
plot(0,0,'*')
xlim([-500 500])
xlabel('x (m)');
ylabel('y (m)');
grid on
title('ray graph')

%print('ray_solution.png','-dpng')
for m=1:length(y0)
    xp(:,m)=[-R-100,x0(m),x2(m),x3(m)];
    yp(:,m)=[y0(m),y0(m),y2(m),y3(m)];
end

ng=zeros(Nx,Ny);
sEfac=zeros(Nx,Ny);
Hrel=zeros(Nx,Ny);
scosE=zeros(Nx,Ny);
ssinE=zeros(Nx,Ny);
ee=zeros(Nx,Ny,nthetabin);
nee=zeros(Nx,Ny,nthetabin);
for j=1:Ny;
    for i=1:Nx;
        xg(i,j)=-R-100+(i-1)*dx;
        yg(i,j)=y0(j);
        for m=1:Ny
            r=hypot(xg(i,j),yg(i,j));
            if r>R
                Ef=1;
            else
                Ef=Efac(m);
            end
            if r>R+100
                Ef=nan;
            end
            [mindist]=get_distance_to_polyline(xg(i,j),yg(i,j),xp(:,m)',yp(:,m)');
            if abs(mindist)<critdist
                ng(i,j)=ng(i,j)+1;
                sEfac(i,j)=sEfac(i,j)+Ef;
                scosE(i,j)=scosE(i,j)+Ef*cosd(theta1(m));
                ssinE(i,j)=ssinE(i,j)+Ef*sind(theta1(m));
                ith=floor((theta1(m)-thetabin0)/thetabin)+1;
                ee(i,j,ith)=ee(i,j,ith)+Ef;
                nee(i,j,ith)=nee(i,j,ith)+1;
            end
        end
        Hrel(i,j)=sqrt(sEfac(i,j)/(2*critdist/dy));
        meandir(i,j)=atan2d(ssinE(i,j),scosE(i,j));
        ee(i,j,:)=ee(i,j,:)/thetabin./nee(i,j,:);
    end
end
meandir(Hrel<thr)=nan;
Hrel(Hrel<thr)=nan;
subplot(122)
% pcolor(xg,yg,ng);
% shading interp
% colorbar
% title('ray density')

% figure;
pcolor(xg,yg,Hrel*.1);
shading interp
axis equal
colorbar
hold on
plot(xpt,ypt,'ro')
title('wave height (m)')
xlim([-500 500])
caxis([0 0.4])
xlabel('x (m)');
grid on
print('reef_focus_wave_height.png','-dpng')

r=hypot(xg,yg);
meandir(r>R)=nan;
figure;
pcolor(xg,yg,270-meandir);
shading interp
axis equal
colorbar
title('wave direction')
ylim([-500 500])
print('reef_focus_wave_dir.png','-dpng')
