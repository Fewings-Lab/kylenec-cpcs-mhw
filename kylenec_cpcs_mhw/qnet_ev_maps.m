% qnet_ev_maps.m

% open the daily ERA5 SST' struct over CPS, SST is in units of
% Kelvin
latlim = [-50 -15];
lonlim = 360 + [-90 -70];

[out,info] = synoptic_anomaly_sfc({'netlw_sfc','netsw_sfc','slhf','sshf'},latlim,lonlim,'01/01/1979','12/21/2020','anomaly');

% dat = read_era5_new('sfc','daily','01/01/1979','12/31/2020',{'sst'},lonlim,latlim,'anomaly',1);
% fields = fieldnames(dat); 
% data = eval(join(['dat.',fields{8}],'.')); % outputs data vars retrieved from read_era5_new into a generically-named matrix
% 
% % filter to time scales of 10 days to 6 months with pl66
% 
% % low-pass filter with a 10-day pl66 filter
% dt = dat.time(2)-dat.time(1); % sampling interval in matdate (number of days)
% Tl = 10; % half-amplitude cutoff period in days to match matdate
% dat_low = NaN(size(data)); % this will be the 10-day low-pass filtered data cube
% for i=1:length(dat.lon) % looping through longitude slices since pl66 only handles 2D matrices
%     dat_low(:,i,:) = (pl66tn(squeeze(data(:,i,:)),dt*24,Tl*24))'; % pl66 assumes time in hours
% end
% 
% 
% % high-pass filter with 6-month pl66 filter
% Th = hours(years(0.5));
% dat_lolo = NaN(size(dat_low));
% for j = 1:length(dat.lon)
%     dat_lolo(:,j,:) = (pl66tn(squeeze(dat_low(:,j,:)),dt*24,Th))';
% end
% 
% dat_band = dat_low-dat_lolo;
% dat_band(:,:,1:2*round(Th/(dt*24))) = NaN;
% dat_band(:,:,end-2*round(Th/(dt*24)):end) = NaN;

% put SST' at times t_ev into a data_cube
dat_sum = out.netlw_sfc+out.netsw_sfc+out.slhf+out.sshf;
save('qnet_data.mat','-v7.3')

%%
load qnet_data.mat
load event_dates.mat
ind = ismember(info.netlw_sfc.time,t_dt);
dat_cube = dat_sum(:,:,ind); 

% convert heat flux to temperature change in deg C/day
cp = 3850; % J/kg/degC
rhow = 1025; % kg/m^3
h = 25; % m
% heat flux is in W/m^2 or J/s/m^2 
%multiplicative factor changes to change per day in numerator and changes from heat flux to temperature change in
% denominator
dat_cube_dT = dat_cube.*(24*3600)./(rhow*cp*h);

% apply land-sea-mask b/c surf heat flux includes over land
em = read_era5_new('other','daily','01/01/1979','12/31/2020',{'land_sea_mask'},lonlim,latlim);
% boo = ismissing(em.land_sea_mask,0);
boo = find(round(em.land_sea_mask));
land_mask = zeros(size(dat_cube_dT,[1 2]));
% land_mask(~boo) = NaN;
land_mask(boo) = NaN;
for z = 1:length(squeeze(dat_cube_dT(1,1,:)))
   dat_cube_dT(:,:,z) = dat_cube_dT(:,:,z)+land_mask;
end

save('qnet_ev_data.mat','dat_cube','info','dat_cube_dT')

%% Coastline data set and coordinate limits around Chile-Peru System:
load qnet_ev_data.mat
latlim = [min(info.sshf.lat) max(info.sshf.lat)];
lonlim = [min(info.sshf.lon) max(info.sshf.lon)]; % works for min/max if lon is in [0 360]
load coastlines
Lat = info.sshf.lat*ones(length(info.sshf.lon),1)';
Lon = ones(length(info.sshf.lat),1)*info.sshf.lon';
clim = [min(dat_cube_dT,[],'all') max(dat_cube_dT,[],'all')]; % limits for the colorbar

%% Plot filled contours of net surface heat flux on world map for each event
load qnet_ev_data.mat
load event_dates.mat

t_dt = rmmissing(t_dt);

for i = 1:ceil(length(t_dt)/16)
   figure()
   for j = 1:16
       k = j+(i-1)*16;
      if k<=length(t_dt)
            subplot(4,4,j)
            h = worldmap(latlim,lonlim); % Map over Chile-Peru System
            setm(h,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
            plotm(coastlat,coastlon) % Adds coastlines
            
            [C,~] = contourm(Lat,Lon,dat_cube_dT(:,:,k),100,'Fill','on'); % Contour of SST'
            caxis(clim)
            cmocean('balance','pivot',0)
            title(datestr(t_dt(k),24),'FontSize',14)
            daspect([1,1.75,1])
            tightmap
            % patchm([-36.5 -35.5 -35.5 -36.5],[-74.25 -74.25 -73.25 -73.25],'k','FaceAlpha',0.5,'LineStyle','--')
      end
      
   end
   if k>length(t_dt)
    sgtitle("Daily Band-pass Filtered Temperature Change from Q'_{net}",'Interpreter','tex','FontSize',20)
    c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
    c.Label.Interpreter = 'latex';
    c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ from $Q'_{net}$ [$^{\circ}$C/day]";
    c.Label.FontSize = 18;
   else
      sgtitle("Daily Band-pass Filtered Q'_{net}",'Interpreter','tex','FontSize',20)
      hp4 = get(subplot(4,4,16),'Position');
      c = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1]);
      c.Label.Interpreter = 'latex';
      c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^{\circ}$C/day]";
      c.Label.FontSize = 18;
   end
end

%% map of conditional average of events
load qnet_ev_data.mat
% load event_dates.mat
load stat_sig_dSST.mat
load 0_05_level_lines.mat

qnet_avg = mean(dat_cube_dT,3,'omitnan');

alpha = 0.05;

N = size(dat_cube_dT,3);

sigma = std(dat_cube_dT,1,3,'omitnan');

p = 1-(alpha/2);
q_t = tinv(p,N-1);

delta_mu = (q_t/sqrt(N)).*sigma;

% mask of statistically significant means
sigmask = xor(qnet_avg>delta_mu,qnet_avg<-delta_mu);

qnet_avg_stat_sig = qnet_avg;
qnet_avg_stat_sig(~sigmask) = NaN;

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
[C_qnet,h_qnet] = contourm(Lat,Lon,qnet_avg_stat_sig,lvls,'Fill','on'); % Contour of SST'
caxis(clim)
cmocean('balance',length(lvls),'pivot',0)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "$$ \frac{ \partial SST'}{ \partial t} $$ [$^\circ$C/day]";
% sgtitle("Average Band-pass Filtered Temperature Change from Q'_{net}",'Interpreter','tex','FontSize',20)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
plotm(y1,x1,'k','LineWidth',0.7)
plotm(y2,x2,'k','LineWidth',0.7)
c.Label.FontSize = 8;
c.FontSize = 8;
set(gcf,'Units','centimeters','Position',[0 0 9.5 11],'PaperUnits','centimeters','PaperPosition',[0 0 9.5 11])
set(gca,'FontSize',8)
setm(gca,'FontSize',8)
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');
PL1.Children.ZData = 4;
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');
PL2.Children.ZData = 5;
t1 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8);

% exportgraphics(gcf,'avg_qnet_stat_sig_v8.png','Resolution',1600)

% saveas(gcf,'avg_qnet_stat_sig_v6.png')