% conditional_avg.m
% Kylene 2022-04-19: not sure why this script exists anymore

load dSSTdt_cube.mat 
load event_dates.mat

dSST_avg = mean(dat_cube,3,'omitnan');

% mapping the average time rate of change of the anomaly
% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(dSST_info.lat) max(dSST_info.lat)];
lonlim = [min(dSST_info.lon) max(dSST_info.lon)]; % works for min/max if lon is in [0 360]
load coastlines
Lat = dSST_info.lat*ones(length(dSST_info.lon),1)';
Lon = ones(length(dSST_info.lat),1)*dSST_info.lon';
clim = [min(dSST_avg,[],'all') max(dSST_avg,[],'all')]; % limits for the colorbar

figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,dSST_avg,100,'Fill','on'); % Contour of SST'
caxis(clim)
cmocean('balance','pivot',0)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^\circ$C/day]";
c.Label.FontSize = 18;
sgtitle("Average Band-pass Filtered $$ \frac{ \partial SST'}{ \partial t} $$",'Interpreter','latex','FontSize',20)

