addpath('..\..')
ncgrid='UGera2cdip_stcroix_net.nc';
samples='mediumv2_UTM.xyz';
snapwave_io_sam2bat(ncgrid,samples)
samples='largev4_UTM.xyz';
snapwave_io_sam2bat(ncgrid,samples)