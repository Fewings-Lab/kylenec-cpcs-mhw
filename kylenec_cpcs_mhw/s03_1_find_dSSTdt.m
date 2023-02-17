% find_dSSTdt.m

latlim = [-50 -15];
lonlim = 360 + [-90 -70];

% open daily dSST/dt anomaly data struct over CPS, SST is in units of
% Kelvin
dat = read_era5_new('sfc','daily','01/01/1979','12/31/2020',{'dsstdt'},lonlim,latlim,'anomaly',1);
fields = fieldnames(dat); 
data = eval(join(['dat.',fields{6}],'.')); % outputs data vars retrieved from read_era5_new into a generically-named matrix

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

% save peak in the dSST'/dt at each event time in a data "cube"
load event_dates.mat
ind = ismember(dat.time,t_dt);
dat_cube = dat_band(:,:,ind).*24.*3600; % converts to deg C/day
dSST_info = rmfield(dat,{fields{5},fields{6}});
dSST_info.time = dat.time(ind);
save('dSSTdt_cube.mat','dat_cube','dSST_info')