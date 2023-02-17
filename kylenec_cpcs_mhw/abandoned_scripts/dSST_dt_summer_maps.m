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
ind(1,:) = find(ismembertol(lat,ptLat)); % columns of data array
ind(2,:) = flip(find(ismembertol(lon,ptLon))); % rows of data array

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
for n = 1:length(sumdSSTdates)
    B = datenum(sumdSSTdates(n));
    tol = datenum(hours(3))/max(abs([A(:);B(:)]));
    
    t0 = ismembertol(A,B,tol); % Finds which 6-hourly time is nearest to summerDates(n)
    if n<=16
        figure(1)
        subplot(4,4,n)
        h = worldmap(latlim,lonlim); % Map over Chile-Peru System
        setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim,'MLabelParallel','south')
        plotm(coastlat,coastlon) % Adds coastlines
        if sum(t0)>1
            [C,~] = contourm(Lat,Lon,mean(dSSTa(:,:,t0),3,'omitnan'),100,'Fill','on'); % Contour of SST'
            caxis(clim)
        else
            [C,~] = contourm(Lat,Lon,dSSTa(:,:,t0),100,'Fill','on'); % Contour of SST'
            caxis(clim)
        end
        title(datestr(sumdSSTdates(n),24),'FontSize',14)
        cmocean('balance','pivot',0)
        scatterm(ptLat,ptLon,20,yc,'filled')
    elseif n>=17
        figure(2)
        subplot(4,4,n-16)
        h = worldmap(latlim,lonlim); % Map over Chile-Peru System
        setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim,'MLabelParallel','south')
        plotm(coastlat,coastlon) % Adds coastlines
        if sum(t0)>1
            [C,~] = contourm(Lat,Lon,mean(dSSTa(:,:,t0),3,'omitnan'),100,'Fill','on'); % Contour of SST'
            caxis(clim)
        else
            [C,~] = contourm(Lat,Lon,dSSTa(:,:,t0),100,'Fill','on'); % Contour of SST'
            caxis(clim)
        end
        title(datestr(sumdSSTdates(n),24),'FontSize',14)
        cmocean('balance','pivot',0)
        scatterm(ptLat,ptLon,20,yc,'filled')
    end
end

figure(1)
sgtitle("Band-pass Filtered $$ \frac{ \partial SST'}{ \partial t} $$",'Interpreter','latex','FontSize',20)
hp4 = get(subplot(4,4,16),'Position');
c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^\circ$C/day]";
c.Label.Interpreter = 'latex';
c.Label.FontSize = 18;

figure(2)
% sgtitle("Band-pass Filtered $$ \frac{ \partial SST'}{ \partial t} $$",'Interpreter','latex','FontSize',20)
hp4 = get(subplot(4,4,15),'Position');
c = colorbar('Position', [hp4(1)+2*hp4(3)+0.07  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^\circ$C/day]";
c.Label.Interpreter = 'latex';
c.Label.FontSize = 18;

%% 10-year chunks of time series of bandpass SST' and dSST'dt 

sstAPL = sstSwA(ind(2,1),ind(1,1),:); % SST' at the yellow dot

dsstPL = dSSTa(ind(2,1),ind(1,1),:); % dSST'dt at the yellow dot

% tenyr = round(hour(years(10))/6); % number of samples in 10 years to limit plots

figure(3)
plot(time1,dsstPL)
ylabel("$$ \frac{ \partial SST'}{ \partial t} $$ [$^\circ$C/day]",'Interpreter','latex')
xlabel("date")