% comparets.m
% Kylene Cooley
% 5 Apr 2021
% as script to comp[are the time series at point lavapie in the different
% bandpass regions

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

dn = datenum(time1); % convert to datenum

[sstSw0,sig0] = clim1y3d(sstSw1, dn, 6, 240); % use 3D data cube climatology function

% Take difference between sstLP and sstSw0 to find SST'
sstSwA = sstLP-sstSw0;

% Obtain indices for lat and lon in data arrays
ptLat = [-35.5 -41];
ptLon = [-72.75 -75];
ind = NaN(2);
ind(1,:) = find(ismembertol(lat,ptLat));
ind(2,:) = flip(find(ismembertol(lon,ptLon)));

ts1 = squeeze(sstSwA(ind(2,1),ind(1,1),:)); % timeseries of unfiltered anomaly at point near Punta Lavapie (-35.5, -72.75)

[low1, high1, low2] = bandpassV3(sstSwA,6,240,0.5);

ts2 = squeeze(high1(ind(2,1),ind(1,1),:)); % timeseries of high-pass filtered anomaly at point near Punta Lavapie (-35.5, -72.75)
ts3 = squeeze(low2(ind(2,1),ind(1,1),:)); % timeseries of low-pass filtered anomaly at point near Punta Lavapie (-35.5, -72.75)

clear sstSw0 sstSw1 sstSwA sstLP % make space in ram for bandpass filtered anomaly

BP = low1-low2;
% Replace one window-length with NaNs on each end
BP(1:2*round(hours(years(0.5))))=NaN;
BP(end-2*round(hours(years(0.5))):end)=NaN;

ts4 = squeeze(BP(ind(2,1),ind(1,1),:)); % timeseries of bandpass filtered anomaly at point near Punta Lavapie (-35.5, -72.75)
