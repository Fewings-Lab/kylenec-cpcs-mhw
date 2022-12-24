% dSSTdt_maps_new.m

load dSSTdt_cube.mat

% mapping the maximum time rate of change of the anomaly at each time found
% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(dSST_info.lat) max(dSST_info.lat)];
lonlim = [min(dSST_info.lon) max(dSST_info.lon)]; % works for min/max if lon is in [0 360]
load coastlines
Lat = dSST_info.lat*ones(length(dSST_info.lon),1)';
Lon = ones(length(dSST_info.lat),1)*dSST_info.lon';
clim = [min(dat_cube,[],'all') max(dat_cube,[],'all')]; % limits for the colorbar

% Plot filled contours of dSST'/dt on world map for each event
for i = 1:ceil(length(dSST_info.time)/16)
    figure()
   for j = 1:16
       k = j+(i-1)*16;
      if k<=length(dSST_info.time)
            subplot(4,4,j)
            h = worldmap(latlim,lonlim); % Map over Chile-Peru System
            setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
            plotm(coastlat,coastlon) % Adds coastlines
            
            [C,~] = contourm(Lat,Lon,dat_cube(:,:,k),100,'Fill','on'); % Contour of SST'
            caxis(clim)
            cmocean('balance','pivot',0)
            title(datestr(dSST_info.time(k),24),'FontSize',14)
            daspect([1,1.75,1])
      end
      
   end
   if k>length(dSST_info.time)
      sgtitle("Daily Band-pass Filtered $$ \frac{ \partial SST'}{ \partial t} $$",'Interpreter','latex','FontSize',20)
      c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
      c.Label.Interpreter = 'latex';
      c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^\circ$C/day]";
      c.Label.FontSize = 18;
   else
      sgtitle("Daily Band-pass Filtered $$ \frac{ \partial SST'}{ \partial t} $$",'Interpreter','latex','FontSize',20)
      hp4 = get(subplot(4,4,16),'Position');
      c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
      c.Label.Interpreter = 'latex';
      c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^\circ$C/day]";
      c.Label.FontSize = 18;
   end
end
