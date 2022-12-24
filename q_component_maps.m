% q_component_maps.m

% parpool(20)

% open the daily ERA5 heat flux data struct saved from Q'net over CPS, in units of
% W/m^2

load qnet_data.mat
load event_dates.mat
clear dat_sum
ind = ismember(info.netlw_sfc.time,t_dt);

netlw_cube = out.netlw_sfc(:,:,ind); 
netsw_cube = out.netsw_sfc(:,:,ind); 
slhf_cube = out.slhf(:,:,ind); 
sshf_cube = out.sshf(:,:,ind); 

% convert heat flux to temperature change in deg C/day
cp = 3850; % J/kg/degC
rhow = 1025; % kg/m^3
h = 25; % m
% heat flux is in W/m^2 or J/s/m^2 
%multiplicative factor changes to change per day in numerator and changes from heat flux to temperature change in
% denominator
netlw_cube_dT = netlw_cube.*(24*3600)./(rhow*cp*h); 
netsw_cube_dT = netsw_cube.*(24*3600)./(rhow*cp*h); 
slhf_cube_dT = slhf_cube.*(24*3600)./(rhow*cp*h); 
sshf_cube_dT = sshf_cube.*(24*3600)./(rhow*cp*h);

% apply land-sea-mask b/c surf heat flux includes over land
em = read_era5_new('other','daily','01/01/1979','12/31/2020',{'land_sea_mask'},lonlim,latlim);
% boo = ismissing(em.land_sea_mask,0);
boo = find(round(em.land_sea_mask));
land_mask = zeros(size(sshf_cube_dT,[1 2]));
% land_mask(~boo) = NaN;
land_mask(boo) = NaN;

for w = 1:length(squeeze(netlw_cube_dT(1,1,:)))
   netlw_cube_dT(:,:,w) = netlw_cube_dT(:,:,w)+land_mask;
end
for x = 1:length(squeeze(netsw_cube_dT(1,1,:)))
   netsw_cube_dT(:,:,x) = netsw_cube_dT(:,:,x)+land_mask;
end
for y = 1:length(squeeze(slhf_cube_dT(1,1,:)))
   slhf_cube_dT(:,:,y) = slhf_cube_dT(:,:,y)+land_mask;
end
for z = 1:length(squeeze(sshf_cube_dT(1,1,:)))
   sshf_cube_dT(:,:,z) = sshf_cube_dT(:,:,z)+land_mask;
end

save('q_component_ev_data.mat','*_cube','info','*_cube_dT','t_dt')

%% Plot filled contours of temperature change from net longwave radiation on world map for each event

% Coastline data set and coordinate limits around Chile-Peru System:
load q_component_ev_data.mat
latlim = [min(info.sshf.lat) max(info.sshf.lat)];
lonlim = [min(info.sshf.lon) max(info.sshf.lon)]; % works for min/max if lon is in [0 360]
load coastlines
Lat = info.sshf.lat*ones(length(info.sshf.lon),1)';
Lon = ones(length(info.sshf.lat),1)*info.sshf.lon';
clim = [min(netlw_cube_dT,[],'all') max(netlw_cube_dT,[],'all')]; % limits for the colorbar

t_dt = rmmissing(t_dt);
% figure()
% subplot(4,4,16)
% hp4 = get(subplot(4,4,16),'Position');
% close

% parpool(20)
for i = 1:ceil(length(t_dt)/16)
   figure()
   for j = 1:16
       k = j+(i-1)*16;
      if k<=length(t_dt)
            subplot(4,4,j)
            h = worldmap(latlim,lonlim); % Map over Chile-Peru System
            setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
            plotm(coastlat,coastlon) % Adds coastlines
            
            [C,~] = contourm(Lat,Lon,netlw_cube_dT(:,:,k),100,'Fill','on'); % Contour of SST'
            caxis(clim)
            cmocean('balance','pivot',0)
            title(datestr(t_dt(k),24),'FontSize',14)
            daspect([1,1.75,1])
            tightmap
            % patchm([-36.5 -35.5 -35.5 -36.5],[-74.25 -74.25 -73.25 -73.25],'k','FaceAlpha',0.5,'LineStyle','--')
      end
      
   end
   if k>length(t_dt)
    sgtitle("Daily Band-pass Filtered Temperature Change from Q'_{LWR}",'Interpreter','tex','FontSize',20)
    c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
    c.Label.Interpreter = 'latex';
    c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^{\circ}$C/day]";
    c.Label.FontSize = 18;
   else
      sgtitle("Daily Band-pass Filtered Temperature Change from Q'_{LWR}",'Interpreter','tex','FontSize',20)
      hp4 = get(subplot(4,4,16),'Position');
      c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
      c.Label.Interpreter = 'latex';
      c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^{\circ}$C/day]";
      c.Label.FontSize = 18;
   end
end

%% map of conditional average of events ofr Q'_LW
load q_component_ev_data.mat
% load event_dates.mat
load stat_sig_dSST.mat

qlw_avg = mean(netlw_cube_dT,3,'omitnan');

alpha = 0.05;

N = size(netlw_cube_dT,3);

sigma = std(netlw_cube_dT,1,3,'omitnan');

p = 1-(alpha/2);
q_t = tinv(p,N-1);

delta_mu = (q_t/sqrt(N)).*sigma;

% mask of statistically significant means
sigmask = xor(qlw_avg>delta_mu,qlw_avg<-delta_mu);

qlw_avg_stat_sig = qlw_avg;
qlw_avg_stat_sig(~sigmask) = NaN;

% mapping the average time rate of change of the anomaly
% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(info.sshf.lat) max(info.sshf.lat)];
lonlim = [min(info.sshf.lon) max(info.sshf.lon)]; % works for min/max if lon is in [0 360]
load coastlines
Lat = info.sshf.lat*ones(length(info.sshf.lon),1)';
Lon = ones(length(info.sshf.lat),1)*info.sshf.lon';
clim = [min(dSST_avg_stat_sig,[],'all') max(dSST_avg_stat_sig,[],'all')]; % limits for the colorbar
lvls = [min(clim,[],'omitnan'):.001:max(clim,[],'omitnan')];

figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,qlw_avg_stat_sig,lvls,'Fill','on'); % Contour of SST'
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^\circ$C/day]";
c.Label.FontSize = 18;
sgtitle("Average Band-pass Filtered Temperature Change from Q'_{LWR}",'Interpreter','tex','FontSize',20)

%% Plot filled contours of temperature change from net shortwave radiation on world map for each event

% Coastline data set and coordinate limits around Chile-Peru System:
load q_component_ev_data.mat
latlim = [min(info.sshf.lat) max(info.sshf.lat)];
lonlim = [min(info.sshf.lon) max(info.sshf.lon)]; % works for min/max if lon is in [0 360]
load coastlines
Lat = info.sshf.lat*ones(length(info.sshf.lon),1)';
Lon = ones(length(info.sshf.lat),1)*info.sshf.lon';
clim = [min(netsw_cube_dT,[],'all') max(netsw_cube_dT,[],'all')]; % limits for the colorbar

t_dt = rmmissing(t_dt);
% figure()
% subplot(4,4,16)
% hp4 = get(subplot(4,4,16),'Position');
% close

% parpool(20)
for i = 1:ceil(length(t_dt)/16)
   figure()
   for j = 1:16
       k = j+(i-1)*16;
      if k<=length(t_dt)
            subplot(4,4,j)
            h = worldmap(latlim,lonlim); % Map over Chile-Peru System
            setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
            plotm(coastlat,coastlon) % Adds coastlines
            
            [C,~] = contourm(Lat,Lon,netsw_cube_dT(:,:,k),100,'Fill','on'); % Contour of SST'
            caxis(clim)
            cmocean('balance','pivot',0)
            title(datestr(t_dt(k),24),'FontSize',14)
            daspect([1,1.75,1])
            tightmap
            % patchm([-36.5 -35.5 -35.5 -36.5],[-74.25 -74.25 -73.25 -73.25],'k','FaceAlpha',0.5,'LineStyle','--')
      end
      
   end
   if k>length(t_dt)
    sgtitle("Daily Band-pass Filtered Temperature Change from Q'_{SWR}",'Interpreter','tex','FontSize',20)
    c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
    c.Label.Interpreter = 'latex';
    c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^{\circ}$C/day]";
    c.Label.FontSize = 18;
   else
      sgtitle("Daily Band-pass Filtered Temperature Change from Q'_{SWR}",'Interpreter','tex','FontSize',20)
      hp4 = get(subplot(4,4,16),'Position');
      c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
      c.Label.Interpreter = 'latex';
      c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^{\circ}$C/day]";
      c.Label.FontSize = 18;
   end
end

%% map of conditional average of events ofr Q'_short wave
load q_component_ev_data.mat
% load event_dates.mat
load stat_sig_dSST.mat

qsw_avg = mean(netsw_cube_dT,3,'omitnan');

alpha = 0.05;

N = size(netsw_cube_dT,3);

sigma = std(netsw_cube_dT,1,3,'omitnan');

p = 1-(alpha/2);
q_t = tinv(p,N-1);

delta_mu = (q_t/sqrt(N)).*sigma;

% mask of statistically significant means
sigmask = xor(qsw_avg>delta_mu,qsw_avg<-delta_mu);

qsw_avg_stat_sig = qsw_avg;
qsw_avg_stat_sig(~sigmask) = NaN;

% mapping the average time rate of change of the anomaly
% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(info.sshf.lat) max(info.sshf.lat)];
lonlim = [min(info.sshf.lon) max(info.sshf.lon)]; % works for min/max if lon is in [0 360]
load coastlines
Lat = info.sshf.lat*ones(length(info.sshf.lon),1)';
Lon = ones(length(info.sshf.lat),1)*info.sshf.lon';
clim = [min(dSST_avg_stat_sig,[],'all') max(dSST_avg_stat_sig,[],'all')]; % limits for the colorbar
lvls = [min(clim):.001:max(clim)];


figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,qsw_avg_stat_sig,lvls,'Fill','on'); % Contour of SST'
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^\circ$C/day]";
c.Label.FontSize = 18;
sgtitle("Average Band-pass Filtered Temperature Change from Q'_{SWR}",'Interpreter','tex','FontSize',20)

%% Plot filled contours of temperature change from suface latent heat flux on world map for each event

% Coastline data set and coordinate limits around Chile-Peru System:
load q_component_ev_data.mat
latlim = [min(info.sshf.lat) max(info.sshf.lat)];
lonlim = [min(info.sshf.lon) max(info.sshf.lon)]; % works for min/max if lon is in [0 360]
load coastlines
Lat = info.sshf.lat*ones(length(info.sshf.lon),1)';
Lon = ones(length(info.sshf.lat),1)*info.sshf.lon';
clim = [min(slhf_cube_dT,[],'all') max(slhf_cube_dT,[],'all')]; % limits for the colorbar

t_dt = rmmissing(t_dt);
% figure()
% subplot(4,4,16)
% hp4 = get(subplot(4,4,16),'Position');
% close

% parpool(20)
for i = 1:ceil(length(t_dt)/16)
   figure()
   for j = 1:16
       k = j+(i-1)*16;
      if k<=length(t_dt)
            subplot(4,4,j)
            h = worldmap(latlim,lonlim); % Map over Chile-Peru System
            setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
            plotm(coastlat,coastlon) % Adds coastlines
            
            [C,~] = contourm(Lat,Lon,slhf_cube_dT(:,:,k),100,'Fill','on'); % Contour of SST'
            caxis(clim)
            cmocean('balance','pivot',0)
            title(datestr(t_dt(k),24),'FontSize',14)
            daspect([1,1.75,1])
            tightmap
            % patchm([-36.5 -35.5 -35.5 -36.5],[-74.25 -74.25 -73.25 -73.25],'k','FaceAlpha',0.5,'LineStyle','--')
      end
      
   end
   if k>length(t_dt)
    sgtitle("Daily Band-pass Filtered Temperature Change from Q'_{LHF}",'Interpreter','tex','FontSize',20)
    c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
    c.Label.Interpreter = 'latex';
    c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^{\circ}$C/day]";
    c.Label.FontSize = 18;
   else
      sgtitle("Daily Band-pass Filtered Temperature Change from Q'_{LHF}",'Interpreter','tex','FontSize',20)
      hp4 = get(subplot(4,4,16),'Position');
      c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
      c.Label.Interpreter = 'latex';
      c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^{\circ}$C/day]";
      c.Label.FontSize = 18;
   end
end

%% map of conditional average of events ofr Q'_long wave
load q_component_ev_data.mat
% load event_dates.mat
load stat_sig_dSST.mat

qslhf_avg = mean(slhf_cube_dT,3,'omitnan');

alpha = 0.05;

N = size(slhf_cube_dT,3);

sigma = std(slhf_cube_dT,1,3,'omitnan');

p = 1-(alpha/2);
q_t = tinv(p,N-1);

delta_mu = (q_t/sqrt(N)).*sigma;

% mask of statistically significant means
sigmask = xor(qslhf_avg>delta_mu,qslhf_avg<-delta_mu);

qslhf_avg_stat_sig = qslhf_avg;
qslhf_avg_stat_sig(~sigmask) = NaN;

% mapping the average time rate of change of the anomaly
% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(info.sshf.lat) max(info.sshf.lat)];
lonlim = [min(info.sshf.lon) max(info.sshf.lon)]; % works for min/max if lon is in [0 360]
load coastlines
Lat = info.sshf.lat*ones(length(info.sshf.lon),1)';
Lon = ones(length(info.sshf.lat),1)*info.sshf.lon';
clim = [min(dSST_avg_stat_sig,[],'all') max(dSST_avg_stat_sig,[],'all')]; % limits for the colorbar
lvls = [min(clim,[],'omitnan'):.001:max(clim,[],'omitnan')];


figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,qslhf_avg_stat_sig,lvls,'Fill','on'); % Contour of SST'
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^\circ$C/day]";
c.Label.FontSize = 18;
sgtitle("Average Band-pass Filtered Temperature Change from Q'_{LHF}",'Interpreter','tex','FontSize',20)

%% Plot filled contours of temperature change from sensible heat flux on world map for each event

% Coastline data set and coordinate limits around Chile-Peru System:
load q_component_ev_data.mat
latlim = [min(info.sshf.lat) max(info.sshf.lat)];
lonlim = [min(info.sshf.lon) max(info.sshf.lon)]; % works for min/max if lon is in [0 360]
load coastlines
Lat = info.sshf.lat*ones(length(info.sshf.lon),1)';
Lon = ones(length(info.sshf.lat),1)*info.sshf.lon';
clim = [min(sshf_cube_dT,[],'all') max(sshf_cube_dT,[],'all')]; % limits for the colorbar

t_dt = rmmissing(t_dt);
% figure()
% subplot(4,4,16)
% hp4 = get(subplot(4,4,16),'Position');
% close

% parpool(20)
for i = 1:ceil(length(t_dt)/16)
   figure()
   for j = 1:16
       k = j+(i-1)*16;
      if k<=length(t_dt)
            subplot(4,4,j)
            h = worldmap(latlim,lonlim); % Map over Chile-Peru System
            setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
            plotm(coastlat,coastlon) % Adds coastlines
            
            [C,~] = contourm(Lat,Lon,sshf_cube_dT(:,:,k),100,'Fill','on'); % Contour of SST'
            caxis(clim)
            cmocean('balance','pivot',0)
            title(datestr(t_dt(k),24),'FontSize',14)
            daspect([1,1.75,1])
            tightmap
            % patchm([-36.5 -35.5 -35.5 -36.5],[-74.25 -74.25 -73.25 -73.25],'k','FaceAlpha',0.5,'LineStyle','--')
      end
      
   end
   if k>length(t_dt)
    sgtitle("Daily Band-pass Filtered Temperature Change from Q'_{SHF}",'Interpreter','tex','FontSize',20)
    c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
    c.Label.Interpreter = 'latex';
    c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^{\circ}$C/day]";
    c.Label.FontSize = 18;
   else
      sgtitle("Daily Band-pass Filtered Temperature Change from Q'_{SHF}",'Interpreter','tex','FontSize',20)
      hp4 = get(subplot(4,4,16),'Position');
      c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
      c.Label.Interpreter = 'latex';
      c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^{\circ}$C/day]";
      c.Label.FontSize = 18;
   end
end

%% map of conditional average of events of Q'_sensible heat flux
load q_component_ev_data.mat
% load event_dates.mat
load stat_sig_dSST.mat

qsshf_avg = mean(sshf_cube_dT,3,'omitnan');

alpha = 0.05;

N = size(sshf_cube_dT,3);

sigma = std(sshf_cube_dT,1,3,'omitnan');

p = 1-(alpha/2);
q_t = tinv(p,N-1);

delta_mu = (q_t/sqrt(N)).*sigma;

% mask of statistically significant means
sigmask = xor(qsshf_avg>delta_mu,qsshf_avg<-delta_mu);

qsshf_avg_stat_sig = qsshf_avg;
qsshf_avg_stat_sig(~sigmask) = NaN;

% mapping the average time rate of change of the anomaly
% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(info.sshf.lat) max(info.sshf.lat)];
lonlim = [min(info.sshf.lon) max(info.sshf.lon)]; % works for min/max if lon is in [0 360]
load coastlines
Lat = info.sshf.lat*ones(length(info.sshf.lon),1)';
Lon = ones(length(info.sshf.lat),1)*info.sshf.lon';
clim = [min(dSST_avg_stat_sig,[],'all') max(dSST_avg_stat_sig,[],'all')]; % limits for the colorbar
lvls = [min(clim,[],'omitnan'):.001:max(clim,[],'omitnan')];


figure()
h = worldmap(latlim,lonlim);
setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,qsshf_avg_stat_sig,lvls,'Fill','on'); % Contour of SST'
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^\circ$C/day]";
c.Label.FontSize = 18;
sgtitle("Average Band-pass Filtered Temperature Change from Q'_{SHF}",'Interpreter','tex','FontSize',20)
