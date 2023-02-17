% s04_0_SST_maps_new.m

% Edited 13 Nov 2022 by KMC

load event_dates.mat

% open the daily ERA5 SST' struct over CPS, SST is in units of
% Kelvin
latlim = [-50 -15];
lonlim = 360 + [-90 -70];

dat = read_era5_new('sfc','daily','01/01/1979','12/31/2020',{'sst'},lonlim,latlim,'anomaly',1);
fields = fieldnames(dat); 
data = eval(join(['dat.',fields{8}],'.')); % outputs data vars retrieved from read_era5_new into a generically-named matrix

% filter to time scales of 10 days to 6 months with pl66

% low-pass filter with a 10-day pl66 filter
dt = dat.time(2)-dat.time(1); % sampling interval in matdate (number of days)
Tl = 10; % half-amplitude cutoff period in days to match matdate
dat_low = NaN(size(data)); % this will be the 10-day low-pass filtered data cube
for i=1:length(dat.lon) % looping through longitude slices since pl66 only handles 2D matrices
    dat_low(:,i,:) = (pl66tn(squeeze(data(:,i,:)),dt*24,Tl*24))'; % pl66 assumes time in hours
end


% high-pass filter with 6-month pl66 filter
Th = hours(years(0.5));
dat_lolo = NaN(size(dat_low));
for j = 1:length(dat.lon)
    dat_lolo(:,j,:) = (pl66tn(squeeze(dat_low(:,j,:)),dt*24,Th))';
end

dat_band = dat_low-dat_lolo;
dat_band(:,:,1:2*round(Th/(dt*24))) = NaN;
dat_band(:,:,end-2*round(Th/(dt*24)):end) = NaN;

% put SST' at times t_ev into a data_cube
ind = ismember(dat.time,t_ev);
dat_cube = dat_band(:,:,ind); % the anomaly should already be in units of deg C

save('SST_ev_data.mat','dat_cube','dat','data')

%% mapping the maximum anomaly at each time found in t_ev
load event_dates.mat
load SST_ev_data.mat

no_dSST = isnan(t_dt);
SST_cube = dat_cube(:,:,~no_dSST);
lat = dat.lat;
lon = dat.lon;

% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(lat) max(lat)];
lonlim = [min(lon) max(lon)]; % works for min/max if lon is in [0 360]
load coastlines
Lat = lat*ones(length(lon),1)';
Lon = ones(length(lat),1)*lon';
clim = [min(dat_cube,[],'all') max(dat_cube,[],'all')]; % limits for the colorbar

%% Plot filled contours of SST' on world map for each event

for i = 1:ceil(length(t_ev)/4)
   figure()
   for j = 1:4
       k = j+(i-1)*4;
      if k<=length(t_ev)
            subplot(2,2,j)
            h = worldmap(latlim,lonlim); % Map over Chile-Peru System
            setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
            plotm(coastlat,coastlon) % Adds coastlines
            
            [C,~] = contourm(Lat,Lon,dat_cube(:,:,k),100,'Fill','on'); % Contour of SST'
            caxis(clim)
            cmocean('balance','pivot',0)
            fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
            title(datestr(t_ev(k),24))
            daspect([1,1.75,1])
            tightmap
            % patchm([-36.5 -35.5 -35.5 -36.5],[-74.25 -74.25 -73.25 -73.25],'k','FaceAlpha',0.5,'LineStyle','--')
      end
      
   end
   if k>length(t_ev)
      % sgtitle("Daily Band-pass Filtered SST'",'Interpreter','latex','FontSize',20)
      c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.03  hp4(2)+hp4(3)*2.1]);
      c.Label.Interpreter = 'latex';
      c.Label.String = "SST' [$^\circ$C]";
      % c.Label.FontSize = 20;
      set(gcf,'PaperPosition',[0 0 6 5])
   else
      % sgtitle("Daily Band-pass Filtered SST'",'Interpreter','latex','FontSize',20)
      hp4 = get(subplot(2,2,4),'Position');
      c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.03  hp4(2)+hp4(3)*2.1]);
      c.Label.Interpreter = 'latex';
      c.Label.String = "SST' [$^\circ$C]";
      % c.Label.FontSize = 20;
      set(gcf,'PaperPosition',[0 0 6 5])
   end
end

%% Composite average map of SST' events

% apply land mask
SST_cube = landmaskcube(SST_cube,lat,lon);

% take avg and get stat sig mask
% average wind stress magnitude anomaly and 95% confidence interval mask
alpha = 0.05;
N = size(SST_cube,3);
[SST_avg_stat_sig,SST_avg] = avg_anomaly_significance_mask(SST_cube,N,alpha);

% map stat sig avg
% titlestr = "Event average 10-day to 6-month bandpass SST anomaly";
titlestr = " ";
clabel = "$\textsf{SST' [}^{\circ}\textsf{C]}$";
cmap = 'balance';
[C_SST,h_SST,c_SST] = avgmap(SST_avg_stat_sig,lat,lon,titlestr,clabel,cmap);

c_SST.Label.FontSize = 8;
c_SST.FontSize = 8;
set(gcf,'Units','centimeters','Position',[0 0 9.5 11],'PaperUnits','centimeters','PaperPosition',[0 0 9.5 11])
set(gca,'FontSize',8)
setm(gca,'FontSize',8)
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');
PL1.Children.ZData = 4;
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');
PL2.Children.ZData = 5;
t1 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8);

% exportgraphics(gcf,'avg_SST_BP_stat_sig_v4.pdf','Resolution',1600)
% saveas(gcf,'avg_SST_BP_stat_sig_v1.png')
