% varfig6_upwell.m (Copy)
% Kylene Cooley
% 31 Mar 2021
% A script to plot the standard deviations of SST, the annual climatology,
% and anomalies for times during the upwelling season
% Using this make standard deviation maps for the rate of change dSST'/dt
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

dSST1 = (sstSw1(:,:,2:end)-sstSw1(:,:,1:end-1))*4; % multiply by 4 for rate of warming per day
sig1 = std(dSST1,0,3,'omitnan'); % standard deviation of dSST/dt
sig1u = std(dSST1(:,:,tupwell),0,3,'omitnan'); % standard deviation of dSST/dt during upwelling times

clear dSST1

dn = datenum(time1); % convert to datenum

[sstSw0,~,~] = clim1y3dV2(sstSw1, dn, 6, 240); % use 3D data cube climatology function

clear time time1 dn

dSST0 = (sstSw0(:,:,2:end)-sstSw0(:,:,1:end-1))*4; % d/dt(climatology) approximation
sig0 = std(dSST0,0,3,'omitnan'); % std of d/dt(climatology)
sig0u = std(dSST0(:,:,tupwell),0,3,'omitnan'); % std during upwelling season of d/dt(climatology)

clear dSST0

SSTa = sstSw1-sstSw0;
dSST2 = (SSTa(:,:,2:end)-SSTa(:,:,1:end-1))*4;
sig2 = std(dSST2,0,3,'omitnan');
sig2u = std(dSST2(:,:,tupwell),0,3,'omitnan'); % std of the anomaly during upwelling season

% Take difference between sstLP and sstSw0 to find SST'
sstSwA = sstLP-sstSw0;

clear sstLP sstSw0 sstSw1 SSTa dSST2

% Bandpass SST' by applying 6-month high-pass filter
[low1, high, low2] = bandpassV3(sstSwA,6,240,0.5); % outputs bandpassed anomaly and std of the low-pass and high-pass parts of the signal respectively
% Take high-pass part of signal
clear sstSwA 
BP = low1-low2;
% Replace one window-length with NaNs on each end
BP(1:2*round(hours(years(0.5))/6))=NaN;
BP(end-2*round(hours(years(0.5))/6):end)=NaN;
clear low1
    
dSST4 = (high(:,:,2:end)-high(:,:,1:end-1))*4;
sig4 = std(dSST4,0,3,'omitnan');
sig4u = std(dSST4(:,:,tupwell),0,3,'omitnan'); % std of the high-pass filtered part of the anomaly
clear dSST4 high

dSST5 = (BP(:,:,2:end)-BP(:,:,1:end-1))*4;
sig5 = std(dSST5,0,3,'omitnan');
sig5u = std(dSST5(:,:,tupwell),0,3,'omitnan'); % std of the band-pass filtered part of the anomaly
clear dSST5 BP

dSST6 = (low2(:,:,2:end)-low2(:,:,1:end-1))*4;
sig6 = std(dSST6,0,3,'omitnan');
sig6u = std(dSST6(:,:,tupwell),0,3,'omitnan'); % std of the low-pass filtered part of the anomaly
clear dSST6 low2

% clear high BP low1 low2

% sig4 = std(dSST4,0,3,'omitnan');
% sig4u = std(dSST4(:,:,tupwell),0,3,'omitnan'); % std of the high-pass filtered part of the anomaly

% sig5 = std(dSST5,0,3,'omitnan');
% sig5u = std(dSST5(:,:,tupwell),0,3,'omitnan'); % std of the band-pass filtered part of the anomaly

% sig6 = std(dSST6,0,3,'omitnan');
% sig6u = std(dSST6(:,:,tupwell),0,3,'omitnan'); % std of the low-pass filtered part of the anomaly

% clear dSST4 dSST5 dSST6

% save standard deviation of different parts to .mat file
save('std_dSSTdt.mat','sig0','sig1','sig2','sig4','sig5','sig6','lat','lon')
save('std_dSSTdt_up.mat','sig0u','sig1u','sig2u','sig4u','sig5u','sig6u','lat','lon')


%% Plotting the variablity of all times
% clear
% load('stdSSTup.mat')


% plot variability of SST and anomalies in 2x3 subplot (all dates)
latlim = [min(lat) max(lat)];
lonlim = [min(lon) max(lon)]; % works for negative longitude, but will have to switch min/max if lon is in [0 360]
load coastlines
Lat = ones(length(lon),1)*lat'; % 81x101 arrays for contour plot
Lon = lon*ones(1,length(lat));
clim1 = [min(sig1,[],'all') 0.5]; % colorbar limits for all except the high pass part
% clim2 = [min(sig5,[],'all') 0.1]; % colorbar limits for SST' slower than 10 days
yc = [255, 255, 0;0, 255, 255]; % color vector for time series pts
ptLat = [-35.5 -41]; % location of time series points
ptLon = [-72.75 -75];

figure(9)

subplot(2,3,1)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig1,100,'Fill','on'); % Contour of std of SST
caxis(clim1)
title('Standard deviation of $$\frac{\partial SST}{\partial t}$$','Interpreter','latex')
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC/day]";
c.Label.FontSize = 14;

subplot(2,3,2)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig0,100,'Fill','on'); % Contour of std of SST
% caxis(clim1)
title('Standard deviation of annual $$ \frac{\partial}{\partial t} $$(SST climatology)','Interpreter','latex')
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC/day]";
c.Label.FontSize = 14;

subplot(2,3,3)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig2,100,'Fill','on'); % Contour of std of SST
caxis(clim1)
title("Standard deviation of $$ \frac{\partial SST'}{\partial t} $$",'Interpreter','latex')
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC/day]";
c.Label.FontSize = 14;

subplot(2,3,4)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig4,100,'Fill','on'); % Contour of std of SST
% caxis(clim)
title("Standard deviation of $$ \frac{\partial SST'}{\partial t} $$ in daily to 10-day band",'Interpreter','latex')
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC/day]";
c.Label.FontSize = 14;

subplot(2,3,5)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig5,100,'Fill','on'); % Contour of std of SST
% caxis(clim2)
title("Standard deviation of $$ \frac{\partial SST'}{\partial t} $$ in 10-day to 6-month band",'Interpreter','latex')
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC/day]";
c.Label.FontSize = 14;

subplot(2,3,6)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig6,100,'Fill','on'); % Contour of std of SST
% caxis(clim2)
title("Standard deviation of $$ \frac{\partial SST'}{\partial t} $$ in 6-month to 40-yr band",'Interpreter','latex')
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC/day]";
c.Label.FontSize = 14;

sgtitle("Variability of $$ \frac{\partial}{\partial t} $$ SST and anomalies during 1980-2019",'Interpreter','latex')

%% Plotting the variablity of upwelling times
% clear
% load('stdSSTup.mat')


% plot variability of SST and anomalies in 2x3 subplot (all dates)
latlim = [min(lat) max(lat)];
lonlim = [min(lon) max(lon)]; % works for negative longitude, but will have to switch min/max if lon is in [0 360]
load coastlines
Lat = ones(length(lon),1)*lat'; % 81x101 arrays for contour plot
Lon = lon*ones(1,length(lat));
clim1 = [min(sig1u,[],'all') 0.5]; % colorbar limits for all except the high pass part
clim2 = [0.3 1]; % colorbar limits for SST' slower than 10 days
yc = [255, 255, 0;0, 255, 255]; % color vector for time series pts
ptLat = [-35.5 -41]; % location of time series points
ptLon = [-72.75 -75];

figure(10)

subplot(2,3,1)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig1u,100,'Fill','on'); % Contour of std of SST
caxis(clim1)
title('Standard deviation of $$\frac{\partial SST}{\partial t}$$','Interpreter','latex')
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC/day]";
c.Label.FontSize = 14;

subplot(2,3,2)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig0u,100,'Fill','on'); % Contour of std of SST
% caxis(clim1)
title('Standard deviation of annual $$ \frac{\partial}{\partial t} $$(SST climatology)','Interpreter','latex')
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC/day]";
c.Label.FontSize = 14;

subplot(2,3,3)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig2u,100,'Fill','on'); % Contour of std of SST
caxis(clim1)
title("Standard deviation of $$ \frac{\partial SST'}{\partial t} $$",'Interpreter','latex')
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC/day]";
c.Label.FontSize = 14;

subplot(2,3,4)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig4u,100,'Fill','on'); % Contour of std of SST
% caxis(clim)
title("Standard deviation of $$ \frac{\partial SST'}{\partial t} $$ in daily to 10-day band",'Interpreter','latex')
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC/day]";
c.Label.FontSize = 14;

subplot(2,3,5)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig5u,100,'Fill','on'); % Contour of std of SST
% caxis(clim2)
title("Standard deviation of $$ \frac{\partial SST'}{\partial t} $$ in 10-day to 6-month band",'Interpreter','latex')
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC/day]";
c.Label.FontSize = 14;

subplot(2,3,6)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig6u,100,'Fill','on'); % Contour of std of SST
% caxis(clim2)
title("Standard deviation of $$ \frac{\partial SST'}{\partial t} $$ in 6-month to 40-yr band",'Interpreter','latex')
cmocean('balance')
scatterm(ptLat,ptLon,20,yc,'filled')
c = colorbar();
c.Label.String = "\sigma [^\circC/day]";
c.Label.FontSize = 14;

sgtitle("Variability of $$ \frac{\partial}{\partial t} $$ SST and anomalies during upwelling seasons1980-2019",'Interpreter','latex')


%% map of just 10-d to 6-mo variability

latlim = [min(lat) max(lat)];
lonlim = [min(lon) max(lon)]; % works for negative longitude, but will have to switch min/max if lon is in [0 360]
load coastlines
Lat = ones(length(lon),1)*lat'; % 81x101 arrays for contour plot
Lon = lon*ones(1,length(lat));
clim = [0.08 0.15]; % manually set color axis limits
yc = [255, 255, 0;0, 255, 255]; % color vector for time series pts
ptLat = [-35.5 -41]; % location of time series points
ptLon = [-72.75 -75];

figure(3)
h = worldmap(latlim,lonlim); % Map over Chile-Peru System
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,sig5u,257,'Fill','on'); % Contour of std of SST
caxis(clim)
title("Standard deviation of $$ \frac{\partial SST'}{\partial t} $$ in 10-day to 6-month band during the upwelling season 1980-2020",'Interpreter','latex')
cmap = cat(1,cmocean('thermal',256),[1 0 0]);
colormap(cmap)
scatterm(ptLat,ptLon,20,yc,'filled')
h = patchm([ptLat(1) ptLat(1)-1 ptLat(1)-1 ptLat(1)],[ptLon(1)-(3/5)  ptLon(1)-(3/5)  ptLon(1)-(8/5) ptLon(1)-(8/5)],'EdgeColor','g','FaceColor','none','LineWidth',2);
c = colorbar();
c.Label.String = "\sigma [^\circC/day]";
c.Label.FontSize = 14;