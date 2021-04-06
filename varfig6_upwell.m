% varfig6_upwell.m
% Kylene Cooley
% 31 Mar 2021
% A script to plot the standard deviations of SST, the annual climatology,
% and anomalies for times during the upwelling season
clear

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
time1 = datetime(time1,'ConvertFrom','datenum'); % convert non-interpolated time to datetimes

% Low-pass filter SST with a 10-day pl66 filter
sstLP = NaN(size(sstSw1)); % This will be the 10 low-pass filtered SST swath
for i=1:length(lon) % Looping through longitude slices because the pl66 filter only handles 2D matrices
    sstLP(i,:,:) = (pl66tn(squeeze(sstSw1(i,:,:)),6,240))'; % Don't forget that for 6-hourly data, dt=6
end

load('tupwell.mat','tupwell')

sig1u = std(sstSw1(:,:,tupwell),0,3,'omitnan'); % standard deviation of SST

dn = datenum(time1); % convert to datenum

[sstSw0,~,sig0u] = clim1y3dV2(sstSw1, dn, 6, 240); % use 3D data cube climatology function
% and output standard deviation during the upwelling season

sig2u = std((sstSw1(:,:,tupwell)-sstSw0(:,:,tupwell)),0,3,'omitnan'); % std of the anomaly during upwelling season

% Take difference between sstLP and sstSw0 to find SST'
sstSwA = sstLP-sstSw0;

% Bandpass SST' by applying 6-month high-pass filter
[low1, high, low2] = bandpassV3(sstSwA,6,240,0.5); % outputs bandpassed anomaly and std of the low-pass and high-pass parts of the signal respectively
% Take high-pass part of signal
clear sstSw0 sstSw1 sstSwA sstLP
BP = low1-low2;
% Replace one window-length with NaNs on each end
BP(1:2*round(hours(years(0.5))))=NaN;
BP(end-2*round(hours(years(0.5))):end)=NaN;
    
sig4u = std(high(:,:,tupwell),0,3,'omitnan'); % std of the high-pass filtered part of the anomaly
sig5u = std(BP(:,:,tupwell),0,3,'omitnan'); % std of the band-pass filtered part of the anomaly
sig6u = std(low2(:,:,tupwell),0,3,'omitnan'); % std of the low-pass filtered part of the anomaly

% save standard deviation of different parts to .mat file
save('stdSSTup.mat','sig0u','sig1u','sig2u','sig4u','sig5u','sig6u','lat','lon')


%% Plotting the variablity
clear
load('stdSSTup.mat')


% plot variability of SST and anomalies in 2x3 subplot (all dates)
latlim = [min(lat) max(lat)];
lonlim = [min(lon) max(lon)]; % works for negative longitude, but will have to switch min/max if lon is in [0 360]
load coastlines
Lat = ones(length(lon),1)*lat'; % 81x101 arrays for contour plot
Lon = lon*ones(1,length(lat));
clim1 = [min(sig1u,[],'all') max(sig1u,[],'all')]; % colorbar limits for all except the high pass part
clim2 = [0.3 1]; % colorbar limits for SST' slower than 10 days
yc = [255, 255, 0;0, 255, 255]; % color vector for time series pts
ptLat = [-35.5 -41]; % location of time series points
ptLon = [-72.75 -75];

figure(9)

subplot(2,3,1)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim)
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig1u,100,'Fill','on'); % Contour of std of SST
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
[C,~] = contourm(Lat,Lon,sig0u,100,'Fill','on'); % Contour of std of SST
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
[C,~] = contourm(Lat,Lon,sig2u,100,'Fill','on'); % Contour of std of SST
caxis(clim1)
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
[C,~] = contourm(Lat,Lon,sig4u,100,'Fill','on'); % Contour of std of SST
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
[C,~] = contourm(Lat,Lon,sig5u,100,'Fill','on'); % Contour of std of SST
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
[C,~] = contourm(Lat,Lon,sig6u,100,'Fill','on'); % Contour of std of SST
caxis(clim2)
title("Standard deviation of SST' in 6-month to 40-yr band")
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC]";
c.Label.FontSize = 14;

sgtitle("Variability of SST and anomalies during upwelling seasons 1980-2019")
