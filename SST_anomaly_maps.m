% SST_anomaly_maps.m
% Kylene Cooley 
% 29 Jan 2021 
% A script to maps the spatial extent of strong warm SST' at times we found during the
% strong upwelling season (Dec-Feb)

% Open data from the .nc file
finfo = ncinfo('Data/ChileCoast-SST-6H.nc');
data = ncreads(finfo.Filename);
sstSw1 = squeeze(data.sst(1:81,:,1,:))-273.15;% Getting rid of most of the landmass and convert to degC
% sstSw1 for 'raw swath SST', or as close to raw as ERA-5 can be
lat = double(squeeze(data.latitude));
lon = double(squeeze(data.longitude));
lon = lon(1:81); % Drop most of the landmass
time = squeeze(data.time);

% Convert hours since Jan 01, 1900 on Gregorian calendar to datenum that
% Matlab can use:
time1 = double(time)/24 + datenum('1900-01-01 00:00:00');
% time2 = interp1(1:6:6*length(time1),time1,1:6*length(time1)); % interpolates 6-hrly time to hourly, but not exactly on the hour
% time2 = datetime(time2,'ConvertFrom','datenum'); % convert to datetimes
% time2 = dateshift(time2, 'start', 'hour', 'nearest'); % shifts datetimes to nearest hour
time1 = datetime(time1,'ConvertFrom','datenum'); % convert non-interpolated time to datetimes

% % Save the relevant variables to a .mat file
% save('ChileCoast_SST_6H.mat','lat','lon','sstSw1','time1','summerDates')

% [Ignore this because it requires too much RAM:]
% Linearly interpolate 6-hrly SST (constant daily values) onto hourly grid
% sstSw2 = interp3(lat,lon,time1,sstSw1,lat,lon,datenum(time2));


% Import variables from ChileCost_SST_6H.mat
% load('ChileCoast_SST_6H.mat')
% Contains: sstSw1 in degC, lat [-45 -20], lon [-90 -70], time1 which is a
% 6-hourly datetime vector, and summerDates which is a datetime vector of
% times we want to plot at the end of this script

% Low-pass filter SST with a 10-day pl66 filter
sstLP = NaN(size(sstSw1)); % This will be the 10 low-pass filtered SST swath
for i=1:length(lon) % Looping through longitude slices because the pl66 filter only handles 2D matrices
    sstLP(i,:,:) = (pl66tn(squeeze(sstSw1(i,:,:)),6,240))'; % Don't forget that for 6-hourly data, dt=6
end

sig1 = std(sstSw1,0,3,'omitnan');
%%
% Make a climatological annual cycle for each lat, lon pair
dn = datenum(time1); % convert to datenum

[sstSw0,sig0] = clim1y3d(sstSw1, dn, 6, 240); % use 3D data cube climatology function

%%
% std of unfiltered anomaly for map
sig2 = std((sstSw1-sstSw0),0,3,'omitnan');

% Take difference between sstLP and sstSw0 to find SST'
sstSwA = sstLP-sstSw0;

% Bandpass SST' by applying 6-month high-pass filter
[sstSwA, sig4, sig6] = bandpass(sstSwA,6,240,0.5); % outputs bandpassed anomaly and std of the low-pass and high-pass parts of the signal respectively
sig5 = std(sstSwA,0,3,'omitnan'); % std of the band-pass filtered part of the anomaly

% save standard deviation of different parts to .mat file
save('stdSST.mat','sig0','sig1','sig2','sig4','sig5','sig6')

%% 
% Find the points corresponding to the time series that we will plot or
% "set our clock" with
ptLat = [-35.5 -41];
ptLon = [-72.75 -75];
ind = NaN(2);
ind(1,:) = find(ismembertol(lat,ptLat));
ind(2,:) = flip(find(ismembertol(lon,ptLon)));

% Extract and plot these time series
% note sstSwA dimension order is lon,lat,time and ind dimension order is
% 1st row lat, 2nd row lon
ts1 = squeeze(sstSwA(ind(2,1),ind(1,1),:)); % point near Punta Lavapie (-35.5, -72.75)
ts2 = squeeze(sstSwA(ind(2,2),ind(1,2),:)); % point further south

% Time-lagged cross correlation time series between these two time series
rho12 = crosscorrTL(ts1,ts2,480);% the last argument, or max lag is largest number of samples between two compared points
% if the max lag in time is 120 days of 6-hourly data (4 samples/day) this
% is a max lag of 240 samples

% plot
% figure(1)
% plot(time1,ts1,'k-',time1,ts2,'r-.')
% legend("SST' at "+num2str(abs(ptLat(1)))+"S, "+num2str(abs(ptLon(1)))+"W", "SST' at "+num2str(abs(ptLat(2)))+"S, "+num2str(abs(ptLon(2)))+"W")
% title("Sea surface temperature anomaly at two points in the CPS")
% xlabel("Date")
% ylabel("SST' [^\circC]")

figure(2)
plot((-120:1/4:120),rho12)
title("Cross-correlation of SST' at (-35.5, -72.75) and (-41, -75)")
ylabel('\rho_{1 2}(\tau)')
xlabel('Time Lag [Days]')
%%
% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(lat) max(lat)];
lonlim = [min(lon) max(lon)]; % works for negative longitude, but will have to switch min/max if lon is in [0 360]
load coastlines
Lat = ones(length(lon),1)*lat';
Lon = lon*ones(1,length(lat));
clim = [min(sstSwA,[],'all') max(sstSwA,[],'all')];
yc = [255, 255, 0;0, 255, 255]; 
load('summerDates.mat')

A = datenum(time1);
% Plot contours of SST' on world map for each date in summerDates
for n = 1:length(summerDates)
    B = datenum(summerDates(n));
    tol = datenum(hours(3))/max(abs([A(:);B(:)]));
    
    t0 = ismembertol(A,B,tol); % Finds which 6-hourly time is nearest to summerDates(n)
    if n<=16
        figure(2)
        subplot(4,4,n)
        h = worldmap(latlim,lonlim); % Map over Chile-Peru System
        setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim)
        plotm(coastlat,coastlon) % Adds coastlines
        if sum(t0)>1
            [C,~] = contourm(Lat,Lon,mean(sstSwA(:,:,t0),3,'omitnan'),100,'Fill','on'); % Contour of SST'
            caxis(clim)
        else
            [C,~] = contourm(Lat,Lon,sstSwA(:,:,t0),100,'Fill','on'); % Contour of SST'
            caxis(clim)
        end
        title(datestr(summerDates(n)))
        cmocean('balance','pivot',0)
        scatterm(ptLat,ptLon,20,yc,'filled')
    elseif n>=17
        figure(3)
        subplot(4,4,n-16)
        h = worldmap(latlim,lonlim); % Map over Chile-Peru System
        setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim)
        plotm(coastlat,coastlon) % Adds coastlines
        if sum(t0)>1
            [C,~] = contourm(Lat,Lon,mean(sstSwA(:,:,t0),3,'omitnan'),100,'Fill','on'); % Contour of SST'
            caxis(clim)
        else
            [C,~] = contourm(Lat,Lon,sstSwA(:,:,t0),100,'Fill','on'); % Contour of SST'
            caxis(clim)
        end
        title(datestr(summerDates(n)))
        cmocean('balance','pivot',0)
        scatterm(ptLat,ptLon,20,yc,'filled')
    end
end

figure(2)
sgtitle("Band-pass Filtered Sea Surface Temperature Anomaly")
hp4 = get(subplot(4,4,16),'Position');
c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
c.Label.String = "SST' [^\circC]";
c.Label.FontSize = 14;
figure(3)
sgtitle("Band-pass Filtered Sea Surface Temperature Anomaly")
hp4 = get(subplot(4,4,16),'Position');
c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
c.Label.String = "SST' [^\circC]";
c.Label.FontSize = 14;

%% 
% % Scatter plot of the two time series against each other (Correlation plot)
% 
% figure(4)
% scatter(ts1,ts2,2,'filled')
% title("Correlation of SST' between northern and southern points")
% xlabel("SST' time series 1 (Punta Lavapie)")
% ylabel("SST' time series 2 (southern point)")
% 
% % Selecting events that occur in Dec-Feb
% dv = datevec(dn);
% monthnum = dv(:, 2); % Grab month number vector
% tupwell = find(monthnum <3 | monthnum==12); % corresponding time indices
% 
% figure(5)
% scatter(ts1(tupwell),ts2(tupwell),2,'filled')
% title("Correlation of SST' between northern and southern points during the upwelling season")
% xlabel("SST' time series 1 (Punta Lavapie)")
% ylabel("SST' time series 2 (southern point)")
% 
% 
% 
% % Scatter plot of the dates we mapped
% D = datenum(summerDates);
% tol = datenum(hours(3))/max(abs([A(:);D(:)]));
% E = ismembertol(A,D,tol);
% 
% figure(6)
% scatter(ts1(E),ts2(E),6,'filled')
% title("Correlation of SST' between northern and southern points during dates of max SST' during the upwelling season")
% xlabel("SST' time series 1 (Punta Lavapie)")
% ylabel("SST' time series 2 (southern point)")




%% 
% plot variability of SST and anomalies in 2x3 subplot (all dates)
latlim = [min(lat) max(lat)];
lonlim = [min(lon) max(lon)]; % works for negative longitude, but will have to switch min/max if lon is in [0 360]
load coastlines
Lat = ones(length(lon),1)*lat'; % 81x101 arrays for contour plot
Lon = lon*ones(1,length(lat));
clim1 = [min(sig1,[],'all') 2.75]; % colorbar limits for SST and climatology
clim2 = [0 max(sig2,[],'all')]; % colorbar limits for SST' slower than 10 days
yc = [255, 255, 0;0, 255, 255]; % color vector for time series pts

figure(8)

subplot(2,3,1)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim)
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig1,100,'Fill','on'); % Contour of std of SST
caxis(clim1)
title('Standard deviation of SST')
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC]";
c.Label.FontSize = 14;

subplot(2,3,2)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim)
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig0,100,'Fill','on'); % Contour of std of SST
caxis(clim1)
title('Standard deviation of annual SST climatology')
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC]";
c.Label.FontSize = 14;

subplot(2,3,3)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim)
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig2,100,'Fill','on'); % Contour of std of SST
caxis(clim2)
title("Standard deviation of SST'")
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC]";
c.Label.FontSize = 14;

subplot(2,3,4)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim)
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig4,100,'Fill','on'); % Contour of std of SST
% caxis(clim)
title("Standard deviation of SST' in daily to 10-day band")
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC]";
c.Label.FontSize = 14;

subplot(2,3,5)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim)
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig5,100,'Fill','on'); % Contour of std of SST
caxis(clim2)
title("Standard deviation of SST' in 10-day to 6-month band")
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC]";
c.Label.FontSize = 14;

subplot(2,3,6)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim)
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig6,100,'Fill','on'); % Contour of std of SST
caxis(clim2)
title("Standard deviation of SST' in 6-month to 40-yr band")
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC]";
c.Label.FontSize = 14;

sgtitle("Variability of SST and anomalies Jun 1979-Dec 2019")
