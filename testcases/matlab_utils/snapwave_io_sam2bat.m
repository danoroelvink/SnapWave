function status = snapwave_io_sam2bat(ncgrid,samples)
%% Simple subroutine to interpolate sampleset on nc grid
%
% Johan Reyns, Deltares, 2023
%
% Input:
%    - ncgrid: netcdf grid name
%    - samples: sample set filename

status = 0;
try
    % read samples
    xyz = load(samples);
    
    % Determine which kind of nc file you have. Currently ugrid, classic nc
    if nc_isvar(ncgrid,'NetNode_z')
        % get coordinates
        x = nc_varget(ncgrid,'NetNode_x');
        y = nc_varget(ncgrid,'NetNode_y');
        z = nc_varget(ncgrid,'NetNode_z');
        F=scatteredInterpolant(xyz(:,1),xyz(:,2));
        znew = F(x,y);
        z(isnan(z))=znew(isnan(z));
        nc_varput(ncgrid,'NetNode_z',z);
    end
    
    if nc_isvar(ncgrid,'mesh2d_node_z')
        % get coordinates
        x = nc_varget(ncgrid,'mesh2d_node_x');
        y = nc_varget(ncgrid,'mesh2d_node_y');
        z = nc_varget(ncgrid,'mesh2d_node_z');
        F=scatteredInterpolant(xyz(:,1),xyz(:,2),xyz(:,3),'linear','none');
        znew = F(x,y);
        z(isnan(z))=znew(isnan(z));
        nc_varput(ncgrid,'mesh2d_node_z',z);
    end
    status = 1;
catch
    error('Could not interpolate data on nc grid.')
end
end