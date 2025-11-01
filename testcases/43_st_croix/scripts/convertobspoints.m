xy=load('ERA5_CDIP_points.txt')
[xutm,yutm]=convertCoordinates(xy(:,1),xy(:,2),'CS1.code',4326,'CS2.name','WGS 84 / UTM zone 20N');
out=[xutm,yutm]
save('ERA5_CDIP_points_utm.txt','out','-ascii')

xy=load('obspoints_stcroix.txt')
[xutm,yutm]=convertCoordinates(xy(:,1),xy(:,2),'CS1.code',4326,'CS2.name','WGS 84 / UTM zone 20N');
out=[xutm,yutm]
save('obspoints_stcroix_utm.txt','out','-ascii')