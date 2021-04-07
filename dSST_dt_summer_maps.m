% dSST_dt_summer_maps.m
% Kylene Cooley
% 6 Apr 2021
% a script to load SST anomaly, time, and lat/lon to map the time rate of
% change

clear
load('sstSwA.mat')
load('summerDates.mat')

% 1st order difference approximation to take partial wrt to time of SST anomaly as difference in SST in degC per 6 hours
dSSTa = 4.*(sstSwA(:,:,2:end)-sstSwA(:,:,1:end-1));

% indices of the lat and lon of the yellow and blue points in the data cube
ptLat = [-35.5 -41];
ptLon = [-72.75 -75];
ind = NaN(2);
ind(1,:) = find(ismembertol(lat,ptLat));
ind(2,:) = flip(find(ismembertol(lon,ptLon)));

% mapping the time rate of change of the anomaly
% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(lat) max(lat)];
lonlim = [min(lon) max(lon)]; % works for negative longitude, but will have to switch min/max if lon is in [0 360]
load coastlines
Lat = ones(length(lon),1)*lat';
Lon = lon*ones(1,length(lat));
clim = [min(dSSTa,[],'all') max(dSSTa,[],'all')];
yc = [255, 255, 0;0, 255, 255]; 

A = datenum(time1(1:end-1));
% Plot contours of SST' on world map for each date in summerDates
for n = 1:length(summerDates)
    B = datenum(summerDates(n));
    tol = datenum(hours(3))/max(abs([A(:);B(:)]));
    
    t0 = ismembertol(A,B,tol); % Finds which 6-hourly time is nearest to summerDates(n)
    if n<=16
        figure(1)
        subplot(4,4,n)
        h = worldmap(latlim,lonlim); % Map over Chile-Peru System
        setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim)
        plotm(coastlat,coastlon) % Adds coastlines
        if sum(t0)>1
            [C,~] = contourm(Lat,Lon,mean(dSSTa(:,:,t0),3,'omitnan'),100,'Fill','on'); % Contour of SST'
            caxis(clim)
        else
            [C,~] = contourm(Lat,Lon,dSSTa(:,:,t0),100,'Fill','on'); % Contour of SST'
            caxis(clim)
        end
        title(datestr(summerDates(n)))
        cmocean('balance','pivot',0)
        scatterm(ptLat,ptLon,20,yc,'filled')
    elseif n>=17
        figure(2)
        subplot(4,4,n-16)
        h = worldmap(latlim,lonlim); % Map over Chile-Peru System
        setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim)
        plotm(coastlat,coastlon) % Adds coastlines
        if sum(t0)>1
            [C,~] = contourm(Lat,Lon,mean(dSSTa(:,:,t0),3,'omitnan'),100,'Fill','on'); % Contour of SST'
            caxis(clim)
        else
            [C,~] = contourm(Lat,Lon,dSSTa(:,:,t0),100,'Fill','on'); % Contour of SST'
            caxis(clim)
        end
        title(datestr(summerDates(n)))
        cmocean('balance','pivot',0)
        scatterm(ptLat,ptLon,20,yc,'filled')
    end
end

figure(1)
sgtitle("Band-pass Filtered Sea Surface Temperature Anomaly")
hp4 = get(subplot(4,4,16),'Position');
c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
c.Label.String = "SST' [^\circC]";
c.Label.FontSize = 14;
figure(2)
sgtitle("Band-pass Filtered Sea Surface Temperature Anomaly")
hp4 = get(subplot(4,4,16),'Position');
c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
c.Label.String = "SST' [^\circC]";
c.Label.FontSize = 14;