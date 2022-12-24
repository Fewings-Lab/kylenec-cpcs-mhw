% windstcrl.m
% Kylene Cooley
% Just mapping wind stress curl outputs from L2 scatterometer data provided
% by Larry
% 10 Mar 2022

%% loading scatterometer wind stress curl outputs with lat/lon grids

% variables in each struct:
% djf_mean_strcrl (240x120), ascat25 (240x200)
% event_daily_anomaly_mean_strcrl (240x120), ascat25 (240x200)
% lat (1x240)
% lon (1x120), ascat25 (1x200)

qscat = load('quikscat_mar22_exp8-output.mat'); % QuikSCAT 11/1/1999-10/31/2009
ascat25 = load('ascata_25km_mar22_exp9-output.mat'); % KNMI ASCAT-A 25km 6/1/2007-5/31/2021
ascat12 = load('ascata_coastal_mar22_exp10-output.mat'); % ASCAT-A 12km (coastal processing) 9/1/2010-8/31/2021

% To check bounds of different lon vectors:
lon1 = qscat.lon;
lon2 = ascat25.lon;
% lon3 = ascat12.lon;
% [lon1(1) lon1(end) lon2(1) lon2(end) lon3(1) lon3(end)]
lat = qscat.lat;

% convert single arrays to double type
qscat.event_daily_anomaly_mean_strcrl = double(qscat.event_daily_anomaly_mean_strcrl);
Qstd = double(qscat.event_daily_anomaly_stddev_strcrl);
ascat25.event_daily_anomaly_mean_strcrl = double(ascat25.event_daily_anomaly_mean_strcrl);
A25std = double(ascat25.event_daily_anomaly_stddev_strcrl(:,81:end));
ascat12.event_daily_anomaly_mean_strcrl = double(ascat12.event_daily_anomaly_mean_strcrl);
A12std = double(ascat12.event_daily_anomaly_stddev_strcrl);

% ascat25 shows 20deg lon farther west (240.125 instead of 260.125)
% They all have the same lon interval: dlon=0.25
% use lon2(81) onwards when combining

load('event_dates.mat')

%% map qscat DJF mean separately

% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(lat) max(lat)];
lonlim = [min(lon1) max(lon1)]; % works for min/max if lon is in [0 360]
load coastlines
% Lat = dSST_info.lat*ones(length(dSST_info.lon),1)';
% Lon = ones(length(dSST_info.lat),1)*dSST_info.lon';
% clim = [min(qscat.djf_mean_strcrl,[],'all') max(qscat.djf_mean_strcrl,[],'all')]; % limits for the colorbar
clim = [-3*10^-7 3*10^-7]; % daily wind stress curl is very small, 
% this makes differences in positive or negative values stand out
lvls = [min(clim,[],'omitnan'):(10^-9):max(clim,[],'omitnan')];

figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat,lon1,qscat.djf_mean_strcrl,lvls,'Fill','on'); % Contour of SST'
% plotm(y1,x1,'k')
% plotm(y2,x2,'k')
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "Wind Stress Curl [N m$^{-3}$]";
c.Label.FontSize = 20;
% sgtitle("Average Band-pass Filtered $$ \frac{ \partial SST'}{ \partial t} $$",'Interpreter','latex','FontSize',20)
set(gcf,'PaperPosition',[0 0 5 6])

% saveas(gcf,'djf_mean_qscat_v2.png')

%% Statistically significant QuikSCAT anomaly field
% I forgot to get the number of events in just the QuikSCAT fields
% I'm thinking potentially the better way to do this is to get separate
% statistically significant fields, then average those together but keeping
% a nan any place that one of them has a nan, like with the Argo MLD
% climatology

%% qscat event mean

% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(lat) max(lat)];
lonlim = [min(lon1) max(lon1)]; % works for min/max if lon is in [0 360]
load coastlines
% Lat = dSST_info.lat*ones(length(dSST_info.lon),1)';
% Lon = ones(length(dSST_info.lat),1)*dSST_info.lon';
clim = [min(qscat.event_daily_anomaly_mean_strcrl,[],'all') max(qscat.event_daily_anomaly_mean_strcrl,[],'all')]; % limits for the colorbar
% clim = [-3*10^-7 3*10^-7]; % daily wind stress curl is very small, 
% this makes differences in positive or negative values stand out
lvls = [min(clim,[],'omitnan'):(10^-8):max(clim,[],'omitnan')];

figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat,lon1,qscat.event_daily_anomaly_mean_strcrl,lvls,'Fill','on'); % Contour of SST'
% plotm(y1,x1,'k')
% plotm(y2,x2,'k')
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "Wind Stress Curl [N m$^{-3}$]";
c.Label.FontSize = 20;
% sgtitle("Average Band-pass Filtered $$ \frac{ \partial SST'}{ \partial t} $$",'Interpreter','latex','FontSize',20)
set(gcf,'PaperPosition',[0 0 5 6])

% saveas(gcf,'event_mean_qscat_v1.png')

%% map ascat25 DJF mean separately

% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(lat) max(lat)];
lonlim = [min(lon2) max(lon2)]; % works for min/max if lon is in [0 360]
load coastlines
% Lat = dSST_info.lat*ones(length(dSST_info.lon),1)';
% Lon = ones(length(dSST_info.lat),1)*dSST_info.lon';
% clim = [min(qscat.djf_mean_strcrl,[],'all') max(qscat.djf_mean_strcrl,[],'all')]; % limits for the colorbar
clim = [-3*10^-7 3*10^-7]; % daily wind stress curl is very small, 
% this makes differences in positive or negative values stand out
lvls = [min(clim,[],'omitnan'):(10^-9):max(clim,[],'omitnan')];

figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat,lon2,ascat25.djf_mean_strcrl,lvls,'Fill','on'); % Contour of SST'
% plotm(y1,x1,'k')
% plotm(y2,x2,'k')
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "Wind Stress Curl [N m$^{-3}$]";
c.Label.FontSize = 20;
% sgtitle("Average Band-pass Filtered $$ \frac{ \partial SST'}{ \partial t} $$",'Interpreter','latex','FontSize',20)
set(gcf,'PaperPosition',[0 0 5 6])

% saveas(gcf,'djf_mean_ascat25_v1.png')

%% ascat25 event mean

% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(lat) max(lat)];
lonlim = [min(lon2) max(lon2)]; % works for min/max if lon is in [0 360]
load coastlines
% Lat = dSST_info.lat*ones(length(dSST_info.lon),1)';
% Lon = ones(length(dSST_info.lat),1)*dSST_info.lon';
clim = [min(ascat25.event_daily_anomaly_mean_strcrl,[],'all') max(ascat25.event_daily_anomaly_mean_strcrl,[],'all')]; % limits for the colorbar
% clim = [-3*10^-7 3*10^-7]; % daily wind stress curl is very small, 
% this makes differences in positive or negative values stand out
lvls = [min(clim,[],'omitnan'):(10^-8):max(clim,[],'omitnan')];

figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat,lon2,ascat25.event_daily_anomaly_mean_strcrl,lvls,'Fill','on'); % Contour of SST'
% plotm(y1,x1,'k')
% plotm(y2,x2,'k')
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "Wind Stress Curl [N m$^{-3}$]";
c.Label.FontSize = 20;
% sgtitle("Average Band-pass Filtered $$ \frac{ \partial SST'}{ \partial t} $$",'Interpreter','latex','FontSize',20)
set(gcf,'PaperPosition',[0 0 5 6])

% saveas(gcf,'event_mean_ascat25_v1.png')

%% map ascat12 DJF mean separately

% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(lat) max(lat)];
lonlim = [min(lon1) max(lon1)]; % works for min/max if lon is in [0 360]
load coastlines
% Lat = dSST_info.lat*ones(length(dSST_info.lon),1)';
% Lon = ones(length(dSST_info.lat),1)*dSST_info.lon';
% clim = [min(qscat.djf_mean_strcrl,[],'all') max(qscat.djf_mean_strcrl,[],'all')]; % limits for the colorbar
clim = [-3*10^-7 3*10^-7]; % daily wind stress curl is very small, 
% this makes differences in positive or negative values stand out
lvls = [min(clim,[],'omitnan'):(10^-9):max(clim,[],'omitnan')];

figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat,lon1,ascat12.djf_mean_strcrl,lvls,'Fill','on'); % Contour of SST'
% plotm(y1,x1,'k')
% plotm(y2,x2,'k')
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "Wind Stress Curl [N m$^{-3}$]";
c.Label.FontSize = 20;
% sgtitle("Average Band-pass Filtered $$ \frac{ \partial SST'}{ \partial t} $$",'Interpreter','latex','FontSize',20)
set(gcf,'PaperPosition',[0 0 5 6])

% saveas(gcf,'djf_mean_ascat12_v1.png')

%% ascat12 event mean

% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(lat) max(lat)];
lonlim = [min(lon1) max(lon1)]; % works for min/max if lon is in [0 360]
load coastlines
clim = [min(ascat12.event_daily_anomaly_mean_strcrl,[],'all') max(ascat12.event_daily_anomaly_mean_strcrl,[],'all')]; % limits for the colorbar
% clim = [-3*10^-7 3*10^-7]; % daily wind stress curl is very small, 
% this makes differences in positive or negative values stand out
lvls = [min(clim,[],'omitnan'):(10^-8):max(clim,[],'omitnan')];

figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat,lon1,ascat12.event_daily_anomaly_mean_strcrl,lvls,'Fill','on'); % Contour of SST'
% plotm(y1,x1,'k')
% plotm(y2,x2,'k')
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "Wind Stress Curl [N m$^{-3}$]";
c.Label.FontSize = 20;
% sgtitle("Average Band-pass Filtered $$ \frac{ \partial SST'}{ \partial t} $$",'Interpreter','latex','FontSize',20)
set(gcf,'PaperPosition',[0 0 5 6])

% saveas(gcf,'event_mean_ascat12_v1.png')

%% average across missions and plot djf

djf = (qscat.djf_mean_strcrl+ascat12.djf_mean_strcrl+ascat25.djf_mean_strcrl(:,81:end))./3;

% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(lat(41:180)) max(lat(41:180))];
lonlim = [min(lon1(41:end)) max(lon1(41:end))]; % works for min/max if lon is in [0 360]
load coastlines
% clim = [min(djf,[],'all') max(djf,[],'all')]; % limits for the colorbar
clim = [-3*10^-7 3*10^-7]; % daily wind stress curl is very small, 
% this makes differences in positive or negative values stand out
lvls = [min(clim,[],'omitnan'):(10^-8):max(clim,[],'omitnan')];

figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),djf(41:180,41:end),lvls,'Fill','on'); % Contour of SST'
% plotm(y1,x1,'k')
% plotm(y2,x2,'k')
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "Wind Stress Curl [N m$^{-3}$]";
c.Label.FontSize = 20;
% sgtitle("Average Band-pass Filtered $$ \frac{ \partial SST'}{ \partial t} $$",'Interpreter','latex','FontSize',20)
set(gcf,'PaperPosition',[0 0 5 6])

% saveas(gcf,'djf_mean_allscat_v3.png')

%% mean of event anomalies across satellite missions

event = (qscat.event_daily_anomaly_mean_strcrl+ascat12.event_daily_anomaly_mean_strcrl+ascat25.event_daily_anomaly_mean_strcrl(:,81:end))./3;

% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(lat(41:180)) max(lat(41:180))];
lonlim = [min(lon1(41:end)) max(lon1(41:end))]; % works for min/max if lon is in [0 360]
load coastlines
% clim = [min(djf,[],'all') max(djf,[],'all')]; % limits for the colorbar
clim = [-3*10^-7 3*10^-7]; % daily wind stress curl is very small, 
% this makes differences in positive or negative values stand out
lvls = [min(clim,[],'omitnan'):(10^-8):max(clim,[],'omitnan')];

figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),event(41:180,41:end),lvls,'Fill','on'); % Contour of SST'
% plotm(y1,x1,'k')
% plotm(y2,x2,'k')
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "Wind Stress Curl [N m$^{-3}$]";
c.Label.FontSize = 20;
% sgtitle("Average Band-pass Filtered $$ \frac{ \partial SST'}{ \partial t} $$",'Interpreter','latex','FontSize',20)
set(gcf,'PaperPosition',[0 0 5 6])

% saveas(gcf,'event_mean_allscat_v3.png')

%% Mask over mean values not statistically significant
% Get number of events from data range for each data set
% There was not a separate 12km output file in the email on 15 Mar 2022
% Got standard deviations for ASCAT-A Coastal and added 20 Mar 2022

% For sigma use, Qstd, A25std, A12std

alpha = 0.05;

Nq = sum(t_dt>=datenum(1999,11,1) & t_dt<=datenum(2009,10,31)); % number of events, N, between the start and end dates of satellite product
Na25 = sum(t_dt>=datenum(2007,6,1) & t_dt<=datenum(2021,5,31));
Na12 = sum(t_dt>=datenum(2010,9,1) & t_dt<=datenum(2021,8,31));

p = 1-(alpha/2);
q_tq = tinv(p,Nq-1);
q_ta25 = tinv(p,Na25-1);
q_ta12 = tinv(p,Na12-1);

delta_muQ = (q_tq/sqrt(Nq)).*Qstd;
delta_muA25 = (q_ta25/sqrt(Na25)).*A25std;
delta_muA12 = (q_ta12/sqrt(Na12)).*A12std;

mask_Q = xor(qscat.event_daily_anomaly_mean_strcrl>delta_muQ,qscat.event_daily_anomaly_mean_strcrl<-delta_muQ);
mask_A25 = xor(ascat25.event_daily_anomaly_mean_strcrl(:,81:end)>delta_muA25,ascat25.event_daily_anomaly_mean_strcrl(:,81:end)<-delta_muA25);
mask_A12 = xor(ascat12.event_daily_anomaly_mean_strcrl>delta_muA12,ascat12.event_daily_anomaly_mean_strcrl<-delta_muA12);

qscat_stat_sig = qscat.event_daily_anomaly_mean_strcrl;
qscat_stat_sig(~mask_Q) = NaN;
ascat25_stat_sig = ascat25.event_daily_anomaly_mean_strcrl(:,81:end);
ascat25_stat_sig(~mask_A25) = NaN;
ascat12_stat_sig = ascat12.event_daily_anomaly_mean_strcrl;
ascat12_stat_sig(~mask_A12) = NaN;

%% statistically significant mean of event anomalies across satellite missions
load 0_05_level_lines.mat

% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(lat(41:180)) max(lat(41:180))];
lonlim = [min(lon1(41:end)) max(lon1(41:end))]; % works for min/max if lon is in [0 360]
load coastlines
% clim = [min(djf,[],'all') max(djf,[],'all')]; % limits for the colorbar
clim = [-3*10^-7 3*10^-7]; % daily wind stress curl is very small, 
% this makes differences in positive or negative values stand out
lvls = [min(clim,[],'omitnan'):(10^-8):max(clim,[],'omitnan')];

figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),event_stat_sig(41:180,41:end),lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r')
plotm(y2,x2,'r')
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "Wind Stress Curl [N m$^{-3}$]";
c.Label.FontSize = 20;
% sgtitle("Average Band-pass Filtered $$ \frac{ \partial SST'}{ \partial t} $$",'Interpreter','latex','FontSize',20)
set(gcf,'PaperPosition',[0 0 5 6])

% saveas(gcf,'event_mean_allscat_v4.png')

%% figure with 9 subplots
% For each mission:
% - typical summer wind stress curl
% - stat sig wind stress curl anomaly
% - difference between these (include where stat sig anomaly = 0 or exclude
%   /show in white)
% same clim bounds but only colorbar on right side

load 0_05_level_lines.mat

% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(lat(41:180)) max(lat(41:180))];
lonlim = [min(lon1(41:end)) max(lon1(41:end))]; % works for min/max if lon is in [0 360]
load coastlines
clim = [-3*10^-7 3*10^-7]; % daily wind stress curl is very small, limits for the colorbar 
% this makes differences in positive or negative values stand out
lvls = [min(clim,[],'omitnan'):(10^-9):max(clim,[],'omitnan')];

diffQ = qscat.djf_mean_strcrl+qscat_stat_sig;
diffA25 = ascat25.djf_mean_strcrl(:,81:end)+ascat25_stat_sig;
diffA12 = ascat12.djf_mean_strcrl+ascat12_stat_sig;

figure()

subplot(3,3,1) % map qscat DJF mean separately
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),qscat.djf_mean_strcrl(41:180,41:end),lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r','LineWidth',2)
plotm(y2,x2,'r','LineWidth',2)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)

subplot(3,3,2) % map qscat DJF mean separately
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),qscat_stat_sig(41:180,41:end),lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r','LineWidth',2)
plotm(y2,x2,'r','LineWidth',2)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)

subplot(3,3,3) % map qscat DJF mean separately
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),diffQ(41:180,41:end),lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r','LineWidth',2)
plotm(y2,x2,'r','LineWidth',2)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)

subplot(3,3,4) % map qscat DJF mean separately
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),ascat25.djf_mean_strcrl(41:180,121:end),lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r','LineWidth',2)
plotm(y2,x2,'r','LineWidth',2)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)

subplot(3,3,5) % map qscat DJF mean separately
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),ascat25_stat_sig(41:180,41:end),lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r','LineWidth',2)
plotm(y2,x2,'r','LineWidth',2)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)

subplot(3,3,6) % map qscat DJF mean separately
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),diffA25(41:180,41:end),lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r','LineWidth',2)
plotm(y2,x2,'r','LineWidth',2)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)

subplot(3,3,7) % map qscat DJF mean separately
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),ascat12.djf_mean_strcrl(41:180,41:end),lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r','LineWidth',2)
plotm(y2,x2,'r','LineWidth',2)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)

subplot(3,3,8) % map qscat DJF mean separately
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),ascat12_stat_sig(41:180,41:end),lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r','LineWidth',2)
plotm(y2,x2,'r','LineWidth',2)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)

subplot(3,3,9) % map qscat DJF mean separately
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),diffA12(41:180,41:end),lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r','LineWidth',2)
plotm(y2,x2,'r','LineWidth',2)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)

hp9 = get(subplot(3,3,9),'Position');
c = colorbar('Position', [hp9(1)+hp9(3)+0.01  hp9(2)  0.01  hp9(2)+hp9(3)*3.1]);
c.Label.Interpreter = 'latex';
c.Label.String = "Wind Stress Curl [N m$^{-3}$]";
c.Label.FontSize = 20;

set(gcf,'PaperPosition',[0 0 8 7])

% saveas(gcf,'djfeventdiff_mean_allscat_v4.png')

%% vertical Ekman pumping velocity
% w_Ek = (wind stress curl)/(density*Coriolis parameter)
% Coriolis parameter = 2*omega * sine(lat)
corinv = (2*7.292*10^-5.*sind(lat(41:180))).^-1;
% w_Ek = (2omega/rho_w)*(event)(lat, lon)*sine(lat)(1, lat)

% Difference between DJF climatologies and anomaly (actually a sum) :/
diffQ = qscat.djf_mean_strcrl+qscat_stat_sig;
diffA25 = ascat25.djf_mean_strcrl(:,81:end)+ascat25_stat_sig;
diffA12 = ascat12.djf_mean_strcrl+ascat12_stat_sig;

% qscat DJF mean 
Qdjf = qscat.djf_mean_strcrl(41:180,41:end);
wE_Qdjf = nan(size(Qdjf));
for i = 1:length(lat(41:180))
    wE_Qdjf(i,:) = Qdjf(i,:).*(corinv(i)/1025)*24*3600;
end

% statistically significant mean qscat anomaly
Qsig = qscat_stat_sig(41:180,41:end);
wE_Qsig = nan(size(Qsig));
for i = 1:length(lat(41:180))
    wE_Qsig(i,:) = Qsig(i,:).*(corinv(i)/1025)*24*3600;
end

% qscat djf+anomaly
Qdif = diffQ(41:180,41:end);
wE_Qdif = nan(size(Qdif));
for i = 1:length(lat(41:180))
    wE_Qdif(i,:) = Qdif(i,:).*(corinv(i)/1025)*24*3600;
end

% ascat25 DJF mean 
A25djf = ascat25.djf_mean_strcrl(41:180,121:end);
wE_A25djf = nan(size(A25djf));
for i = 1:length(lat(41:180))
    wE_A25djf(i,:) = A25djf(i,:).*(corinv(i)/1025)*24*3600;
end

% statistically significant mean ascat25 anomaly
A25sig = ascat25_stat_sig(41:180,41:end);
wE_A25sig = nan(size(A25sig));
for i = 1:length(lat(41:180))
    wE_A25sig(i,:) = A25sig(i,:).*(corinv(i)/1025)*24*3600;
end

% ascat25 djf+anomaly
A25dif = diffA25(41:180,41:end);
wE_A25dif = nan(size(A25dif));
for i = 1:length(lat(41:180))
    wE_A25dif(i,:) = A25dif(i,:).*(corinv(i)/1025)*24*3600;
end

% ascat12 DJF mean 
A12djf = ascat12.djf_mean_strcrl(41:180,41:end);
wE_A12djf = nan(size(A12djf));
for i = 1:length(lat(41:180))
    wE_A12djf(i,:) = A12djf(i,:).*(corinv(i)/1025)*24*3600;
end

% statistically significant mean ascat12 anomaly
A12sig = ascat12_stat_sig(41:180,41:end);
wE_A12sig = nan(size(A12sig));
for i = 1:length(lat(41:180))
    wE_A12sig(i,:) = A12sig(i,:).*(corinv(i)/1025)*24*3600;
end

% ascat12 djf+anomaly
A12dif = diffA12(41:180,41:end);
wE_A12dif = nan(size(A12dif));
for i = 1:length(lat(41:180))
    wE_A12dif(i,:) = A12dif(i,:).*(corinv(i)/1025)*24*3600;
end


%% Vertical Ekman pumping velocity anomaly for each scatterometer product 
% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0], [0.1 0.01], [0.01 0.11]);
% if ~make_it_tight,  clear subplot;  end

load 0_05_level_lines.mat
latlim = [min(lat(41:180)) max(lat(41:180))];
lonlim = [min(lon1(41:end)) max(lon1(41:end))]; % works for min/max if lon is in [0 360]
load coastlines
clim = [-5*10^-1 5*10^-1]; % daily vertical Ekman pumping velocity is very small, limits for the colorbar 
% this makes differences in positive or negative values stand out
lvls = [min(clim,[],'omitnan'):(10^-2):max(clim,[],'omitnan')];

figure()
tiledlayout(3,3,'TileSpacing','compact','Padding','compact')

nexttile
%subplot(3,3,1) % map qscat DJF mean 
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','north')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),wE_Qdjf,lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r','LineWidth',0.7)
plotm(y2,x2,'r','LineWidth',0.7)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
setm(h,'FontSize',8)
set(gca,'FontSize',8)
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');
PL1.Children.ZData = 4;
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');
PL2.Children.ZData = 5;
% t1 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8);

nexttile
%subplot(3,3,2) % map statistically significant mean qscat anomaly
h = worldmap(latlim,lonlim);
% setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','north')
setm(h,'ParallelLabel','off','MeridianLabel','off')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),wE_Qsig,lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r','LineWidth',0.7)
plotm(y2,x2,'r','LineWidth',0.7)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
setm(h,'FontSize',8)
set(gca,'FontSize',8)
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');
PL1.Children.ZData = 4;
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');
PL2.Children.ZData = 5;
% t1 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8);

nexttile
%subplot(3,3,3) % map qscat djf+anomaly
h = worldmap(latlim,lonlim);
% setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','north')
setm(h,'ParallelLabel','off','MeridianLabel','off')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),wE_Qdif,lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r','LineWidth',0.7)
plotm(y2,x2,'r','LineWidth',0.7)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
setm(h,'FontSize',8)
set(gca,'FontSize',8)
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');
PL1.Children.ZData = 4;
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');
PL2.Children.ZData = 5;
% t1 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8);

nexttile
%subplot(3,3,4) % map ascat 25 km DJF mean
h = worldmap(latlim,lonlim);
% setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','north')
setm(h,'ParallelLabel','off','MeridianLabel','off')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),wE_A25djf,lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r','LineWidth',0.7)
plotm(y2,x2,'r','LineWidth',0.7)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
setm(h,'FontSize',8)
set(gca,'FontSize',8)
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');
PL1.Children.ZData = 4;
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');
PL2.Children.ZData = 5;
% t1 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8);

nexttile
%subplot(3,3,5) % map statistically significant mean ascat 25 km anomaly
h = worldmap(latlim,lonlim);
% setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','north')
setm(h,'ParallelLabel','off','MeridianLabel','off')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),wE_A25sig,lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r','LineWidth',0.7)
plotm(y2,x2,'r','LineWidth',0.7)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
setm(h,'FontSize',8)
set(gca,'FontSize',8)
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');
PL1.Children.ZData = 4;
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');
PL2.Children.ZData = 5;
% t1 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8);

nexttile
%subplot(3,3,6) % map ascat 25 km djf+anomaly
h = worldmap(latlim,lonlim);
% setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','north')
setm(h,'ParallelLabel','off','MeridianLabel','off')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),wE_A25dif,lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r','LineWidth',0.7)
plotm(y2,x2,'r','LineWidth',0.7)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
setm(h,'FontSize',8)
set(gca,'FontSize',8)
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');
PL1.Children.ZData = 4;
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');
PL2.Children.ZData = 5;
% t1 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8);

nexttile
%subplot(3,3,7) % map ascat 12.5 km DJF mean
h = worldmap(latlim,lonlim);
% setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','north')
setm(h,'ParallelLabel','off','MeridianLabel','off')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),wE_A12djf,lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r','LineWidth',0.7)
plotm(y2,x2,'r','LineWidth',0.7)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
setm(h,'FontSize',8)
set(gca,'FontSize',8)
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');
PL1.Children.ZData = 4;
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');
PL2.Children.ZData = 5;
% t1 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8);

nexttile
%subplot(3,3,8) % map statistically significant mean ascat 12.5 km anomaly
h = worldmap(latlim,lonlim);
% setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','north')
setm(h,'ParallelLabel','off','MeridianLabel','off')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),wE_A12sig,lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r','LineWidth',0.7)
plotm(y2,x2,'r','LineWidth',0.7)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
setm(h,'FontSize',8)
set(gca,'FontSize',8)
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');
PL1.Children.ZData = 4;
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');
PL2.Children.ZData = 5;
% t1 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8);

nexttile
%subplot(3,3,9) % map ascat 12.5 km djf+anomaly
h = worldmap(latlim,lonlim);
% setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','north')
setm(h,'ParallelLabel','off','MeridianLabel','off')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),wE_A12dif,lvls,'Fill','on'); % Contour of SST'
plotm(y1,x1,'r','LineWidth',0.7)
plotm(y2,x2,'r','LineWidth',0.7)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
setm(h,'FontSize',8)
set(gca,'FontSize',8)
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');
PL1.Children.ZData = 4;
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');
PL2.Children.ZData = 5;
% t1 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8);

%hp9 = get(subplot(3,3,9),'Position');
% c = colorbar('Position', [hp9(1)+hp9(3)  0.05  0.01  hp9(2)+hp9(3)*2.7]);
c = colorbar;
c.Layout.Tile = 'east';
c.Label.Interpreter = 'latex';
c.Label.String = "$\textsf{vertical Ekman velocity or anomaly [m day}^{-1}\textsf{]}$";
c.Label.FontSize = 8;
c.FontSize = 8;

set(gcf,'Units','centimeters','Position',[0 0 17 18],'PaperUnits','centimeters','PaperPosition',[0 0 17 18])
% exportgraphics(gcf,'djfeventdiff_mean_wEk_v6.pdf','Resolution',2400)
%set(gcf,'PaperPosition',[0 0 4 5])

%saveas(gcf,'djfeventdiff_mean_wEk_v2.png')


%% map vertical Ekman pumping velocity anomaly

% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(lat(41:180)) max(lat(41:180))];
lonlim = [min(lon1(41:end)) max(lon1(41:end))]; % works for min/max if lon is in [0 360]
load coastlines
% clim = [min(wE,[],'all') max(wE,[],'all')]; % limits for the colorbar
clim = [-5*10^-6 5*10^-6]; % daily wind stress curl is very small, 
% this makes differences in positive or negative values stand out
lvls = [min(clim,[],'omitnan'):(10^-8):max(clim,[],'omitnan')];

figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C_dsst,h_dsst] = contourm(lat(41:180),lon1(41:end),wE(41:180,41:end),lvls,'Fill','on'); % Contour of SST'
% plotm(y1,x1,'k')
% plotm(y2,x2,'k')
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
% plotm(-24,270,'pb','MarkerSize',10,'MarkerFaceColor','b') % to add a small blue star at center of dSST'dt
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
c = colorbar();
c.Label.Interpreter = 'tex';
c.Label.String = "w_{\fontsize{15}\rmEk} [m s^{-1}]";
c.Label.FontSize = 20;
% sgtitle("Average Band-pass Filtered $$ \frac{ \partial SST'}{ \partial t} $$",'Interpreter','latex','FontSize',20)
set(gcf,'PaperPosition',[0 0 5 6])

% saveas(gcf,'wEk_mean_v4.png')