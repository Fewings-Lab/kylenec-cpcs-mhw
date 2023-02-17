% s07_0_MLD_view.m

load Argo_mixedlayers_monthlyclim_03172021.mat
load('windstress_era5_v2.mat','latlim','lonlim')
load coastlines

% Extract the monthly clim for month numbers 1, 2, 12
mld_Dec = squeeze(mld_ta_mean(12,90:110,40:75))';
mld_Jan = squeeze(mld_ta_mean(1,90:110,40:75))';
mld_Feb = squeeze(mld_ta_mean(2,90:110,40:75))';
% average months 1, 2, 12
month_cube = mld_ta_mean([1 2 12],90:110,40:75);
summer_mean = squeeze(mean(month_cube,1,'includenan'))';

% clim = [min(summer_mean,[],'all') max(summer_mean,[],'all')]; % limits for the colorbar
clim = [0 max(summer_mean,[],'all')]; % limits for the colorbar
lvls = [min(clim,[],'omitnan'):1:max(clim,[],'omitnan')];

% plot of Dec climatology
figure()

h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
pcolorm(latm(1,40:75),lonm(90:110,1),mld_Dec)
caxis(clim)
c = colorbar();
cmocean('deep')
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
ylabel('Latitude [^\circN]')
xlabel('Longitude [^\circE]')
c.Label.String = 'Mixed layer depth, h [m]';
c.Label.FontSize = 16;
title('Monthly climatology of MLD in December')


% image(lonm(90:110,1),latm(1,40:75),mld_Dec,'CDataMapping','scaled')
% axis xy
% caxis(clim)
% c = colorbar();
% cmocean('thermal')
% ylabel('Latitude [^\circN]')
% xlabel('Longitude [^\circE]')
% c.Label.String = 'Mixed layer depth, h [m]';
% c.Label.FontSize = 16;
% title('Monthly climatology of MLD in December')

% saveas(gcf,'mld_ta_mean_monthlyclimatology_Dec_v4.png')


% Plot of Jan climatology
figure()

h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
pcolorm(latm(1,40:75),lonm(90:110,1),mld_Jan)
caxis(clim)
c = colorbar();
cmocean('deep')
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
ylabel('Latitude [^\circN]')
xlabel('Longitude [^\circE]')
c.Label.String = 'Mixed layer depth, h [m]';
c.Label.FontSize = 16;
title('Monthly climatology of MLD in January')

% image(lonm(90:110,1),latm(1,40:75),mld_Jan,'CDataMapping','scaled')
% axis xy
% caxis(clim)
% c = colorbar();
% cmocean('thermal')
% ylabel('Latitude [^\circN]')
% xlabel('Longitude [^\circE]')
% c.Label.String = 'Mixed layer depth, h [m]';
% c.Label.FontSize = 16;
% title('Monthly climatology of MLD in January')

% saveas(gcf,'mld_ta_mean_monthlyclimatology_Jan_v4.png')

% Plot of Feb climatology
figure()

h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
pcolorm(latm(1,40:75),lonm(90:110,1),mld_Feb)
caxis(clim)
c = colorbar();
cmocean('deep')
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
ylabel('Latitude [^\circN]')
xlabel('Longitude [^\circE]')
c.Label.String = 'Mixed layer depth, h [m]';
c.Label.FontSize = 16;
title('Monthly climatology of MLD in February')

% image(lonm(90:110,1),latm(1,40:75),mld_Feb,'CDataMapping','scaled')
% axis xy
% caxis(clim)
% c = colorbar();
% cmocean('thermal')
% ylabel('Latitude [^\circN]')
% xlabel('Longitude [^\circE]')
% c.Label.String = 'Mixed layer depth, h [m]';
% c.Label.FontSize = 16;
% title('Monthly climatology of MLD in February')

% saveas(gcf,'mld_ta_mean_monthlyclimatology_Feb_v4.png')

% Plot of summer seasonal MLD climatology

latlim = [min(latm(1,40:75)) max(latm(1,40:75))];
lonlim = [min(lonm(90:110,1)) max(lonm(90:110,1))]; % works for min/max if lon is in [0 360]
load coastlines
% Lat = latm(1,40:75)'*ones(length(lonm(90:110,1)),1)';
% Lon = ones(length(latm(1,40:75)),1)'*lonm(90:110,1)';
load('0_05_level_lines.mat')

figure()
% image(lonm(90:110,1),latm(1,40:75),summer_mean,'CDataMapping','scaled')
% axis xy

h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
pcolorm(latm(1,40:75),lonm(90:110,1),summer_mean)
plotm(y1,x1,'r','LineWidth',2)
plotm(y2,x2,'r','LineWidth',2)
caxis(clim)
c = colorbar();
cmocean('deep')
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
ylabel('Latitude [^\circN]')
xlabel('Longitude [^\circE]')
c.Label.String = 'Mixed layer depth, h [m]';
c.Label.FontSize = 16;
% title('Seasonal climatology of MLD in summer (DJF)')

% saveas(gcf,'mld_ta_mean_monthlyclimatology_summer_v5.png')

% Find a mean value: excluding NaNs and weighting by latitude
lat = latm(1,40:75);
lon = lonm(90:110,1);

% w = sum(cosd(lat));
% X = sum(summer_mean,2,'omitnan')./sum(isfinite(summer_mean),2,'omitnan');
% W = (1/w)*sum(cosd(lat).*X','omitnan');

mld_mean = squeeze(mean(avg1,'omitnan'));

% save MLD climatologies and overall avg
save('mld_mean_eval.mat','mld_Dec','mld_Jan','mld_Feb','summer_mean','mld_mean')