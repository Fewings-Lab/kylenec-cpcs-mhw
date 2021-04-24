% Cross-correlation of dSST'/dt in bandpass
% Kylene Cooley
% 22 Apr 2021
% Calculates the cross-correlation coefficients with the yellow dot point
% (Punta Lavapie) dSST'/dt, maps this and maps this. 95% critical value is bottom of
% the colorbar

load("sstSwA.mat")

dSSTdt = sstSwA(:,:,2:end)-sstSwA(:,:,1:end-1);

% cross-correlation coefficient array
RhoXY = zeros(length(lon),length(lat));

% find and extract the time series at Punta Lavapie
ptLat = -35.5;
ptLon = -72.75;
ind = NaN(1,2);
ind(1) = find(ismembertol(lat,ptLat));
ind(2) = flip(find(ismembertol(lon,ptLon)));
X = dSSTdt(ind(2),ind(1),:);

for j = 1:length(lon)
    for k = 1:length(lat)
        Y = dSSTdt(j,k,:); % time series to compare
        RhoXY(j,k) = squeeze(crosscorrTL(X,Y,0));
    end
end

%% Plot the cross-correlation coefficient map
% not with lower clim = 95% critical value, need help

latlim = [min(lat) max(lat)];
lonlim = [min(lon) max(lon)]; % works for negative longitude, but will have to switch min/max if lon is in [0 360]
load coastlines
Lat = ones(length(lon),1)*lat'; % 81x101 arrays for contour plot
Lon = lon*ones(1,length(lat));
% clim = [0.95 1]; % manually set color axis limits
yc = [255, 255, 0;0, 255, 255]; % color vector for time series pts
ptLat = [-35.5 -41]; % location of time series points
ptLon = [-72.75 -75];

figure(1)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,RhoXY,257,'Fill','on'); % Contour of std of SST
% caxis(clim)
title("Cross-correlation coefficient of $$ \frac{\partial SST'}{\partial t} $$ in 10-day to 6-month band during the upwelling season 1980-2020",'Interpreter','latex')
cmap = cat(1,cmocean('thermal',256),[1 0 0]);
colormap(cmap)
scatterm(ptLat,ptLon,20,yc,'filled')
h = patchm([ptLat(1) ptLat(1)-1 ptLat(1)-1 ptLat(1)],[ptLon(1)-(3/5)  ptLon(1)-(3/5)  ptLon(1)-(8/5) ptLon(1)-(8/5)],'EdgeColor','g','FaceColor','none','LineWidth',2);
c = colorbar();
c.Label.String = "\sigma [^\circC/day]";
c.Label.FontSize = 14;