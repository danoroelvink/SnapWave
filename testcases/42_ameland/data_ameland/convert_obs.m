obsfil='locations_obs.xyn';
WGS_obsfil='locations_WGS_obs.xyn';
[xRD,yRD,statname] = fm_io_getobsfil(obsfil);
[lon,lat]=convertCoordinates(xRD,yRD,'CS1.code',28992,'CS2.code',4326);
fi=fopen(WGS_obsfil,'w')';
for i=1:length(lat)
    fprintf(fi,'%f %f %s \n',lon(i),lat(i),statname{i});
end
fclose(fi)
