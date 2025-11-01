clear all;close all;
addpath('..\..\matlab_utils')
runid={'shoalref_coarse_neu','shoalref_fine_neu','shoalref_refined_neu'};
col={'.b','.r','.g'};
for irun=1:length(runid);
    g=9.81;
    zr=6;b=0.323;A=1.4;Ab=1;xb=300;Rb=200;
    Lb=200;Tb=20;dx=10;t=0;Lgrid=2000;
    dy=10;Lgridy=10000;
    hmin=1;
    z0=0;
    T=5;
    H0=1;
    figure(1);
    subplot(311);
    x=Lgrid:-dx:0;
    zbmean=zr-A*x.^b;
    zb=zbmean-Ab*exp(-((x-xb)/Rb).^2).*cos(2*pi*(x/Lb-t/Tb));
    %h=max(z0-zb,hmin);
    h=z0-zb;
    h(h<=hmin)=nan;
    plot(x,h,'k','linewidth',1);
    axis([150 2000 0 10]);
    ylabel('depth (m)')
    title('Water depth')
    idir=0
    for dir0=[0,45,60]
        idir=idir+1;
        t=0;
        S=zeros(size(x));
        hold on
        [k,C,Cg]=disper(h,T);
        dirm(idir,:)=asind(C/C(1)*sind(dir0));
        Ksh=sqrt(Cg(1)./Cg);
        Kr=sqrt(cosd(dir0)./cosd(dirm(idir,:)));
        Hm(idir,:)=H0.*Ksh.*Kr;
        subplot(312)
        plot(x,dirm(idir,:),'k','linewidth',1)
        xlim([150 2000])
        hold on
        ylabel('dir.(^o)')
        title('Wave direction')
        subplot(313)
        plot(x,Hm(idir,:),'k','linewidth',1)
        xlim([150 2000])
        hold on
        title('Wave height')
        ylabel('H (m)')
        xlabel('cross-shore distance (m)')
    end
    subplot(313)

    nx=Lgrid/dx;
    ny=Lgridy/dy;
    for j=1:ny+1;
        xx(:,j)=x;
        for i=1:nx+1;
            yy(i,j)=(j-1)*dy;
            zz(i,j)=zb(i);
        end
    end
    figure;
    pcolor(xx,yy,zz);
    shading flat;
    colorbar
    out=[xx(:),yy(:),zz(:)];
    save('profile.xyz','out','-ascii');
    obs=[x',x'*0+5000];
    save('central.obs','obs','-ascii');
    fname=['..\output\',runid{irun},'_his.nc'];
    xhis=ncread(fname,'station_x')
    Hm0(irun,:,:)=ncread(fname,'point_hm0');
    wavdir(irun,:,:)=ncread(fname,'point_wavdir');
    figure(1);
    subplot(312)
    for idir=1:3
        plot(xhis,90-squeeze(wavdir(irun,:,idir)),col{irun},'linewidth',1)
    end
    subplot(313)
    for idir=1:3
        plot(xhis,squeeze(Hm0(irun,:,idir)),col{irun},'linewidth',1)
    end
    ylim([0.75 1.25])
    %legend('0^o an','45^o an','60^o an','0^o num','45^o num','60^o num','orientation','vertical')
end
print('..\results\shoalref.png','-r300','-dpng')
for irun=1:length(runid)
    for idir=1:3
        measx=x;measy=squeeze(Hm(idir,:));compx=x;compy=squeeze(Hm0(irun,:,idir));
        [sH(irun,idir)]=skill0([measx(:),measy(:)],[compx(:),compy(:)],'nlast',15);
        measx=x;measy=squeeze(dirm(idir,:));compx=x;compy=squeeze(90-wavdir(irun,:,idir));
        [sdir(irun,idir)]=skill0([measx(:),measy(:)],[compx(:),compy(:)],'nlast',15);
    end
end
dir0=[0,30,45];
sp=sH;
fi=fopen('..\results\shoalref_stats_Hm0.txt','w');
fprintf(fi,'%20s %4s %8s %8s %8s %8s %8s %8s\n',...
    'runid','dir','r2','sci','relbias','bias','skill','rmse')
for idir=1:3
    for irun=1:length(runid)
        fprintf(fi,'%20s %4.0f %8.4f %8.3f %8.3f %8.3f %8.3f %8.3f\n', ...
            runid{irun},dir0(idir), ...
            sp(irun,idir).r2,sp(irun,idir).sci,sp(irun,idir).relbias, ...
            sp(irun,idir).bias,sp(irun,idir).skill,sp(irun,idir).rmse)
    end
end
fclose(fi)

sp=sdir;
fi=fopen('..\results\shoalref_stats_dir.txt','w');
fprintf(fi,'%20s %4s %8s %8s %8s %8s %8s %8s\n',...
    'runid','dir','r2','sci','relbias','bias','skill','rmse')
for idir=2:3
    for irun=1:length(runid)
        fprintf(fi,'%20s %4.0f %8.4f %8.3f %8.3f %8.3f %8.3f %8.3f\n', ...
            runid{irun},dir0(idir), ...
            sp(irun,idir).r2,sp(irun,idir).sci,sp(irun,idir).relbias, ...
            sp(irun,idir).bias,sp(irun,idir).skill,sp(irun,idir).rmse)
    end
end
fclose(fi)

