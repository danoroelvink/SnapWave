function [xmap_out,ymap_out,zmap] = data_density(x,y,width, height)
% Simple function to create data density map
% v1.0  Nederhoff   Jun-18

xmap = linspace(min(x), max(x), width+1);
ymap = linspace(min(y), max(y), height+1);
for ii = 1:length(xmap)-1
    xmap_out(ii) = nanmean(xmap([ii:ii+1]));
    for jj = 1:length(ymap)-1
        if ii ==1;
            ymap_out(jj) = nanmean(ymap([jj:jj+1]));
        end
        zmap(ii,jj) = sum(x>xmap(ii) & x<xmap(ii+1) & y>ymap(jj) & y<ymap(jj+1));
    end
end
end

