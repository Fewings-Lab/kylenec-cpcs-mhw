% SST_anomaly_maps.m
% Kylene Cooley 
% 29 Jan 2021 
% A script to maps the spatial extent of strong warm SST' at times we found during the
% strong upwelling season (Dec-Feb)

% Open data from the .nc file
finfo = ncinfo('Data/ChileCoast-SST-6H.nc');
data = ncreads(finfo.Filename);
sstSw1 = squeeze(data.sst(1:81,:,1,:))-273.15;% Getting rid of most of the landmass
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

% Linearly interpolate 6-hrly SST (constant daily values) onto hourly grid
% sstSw2 = interp3(lat,lon,time1,sstSw1,lat,lon,datenum(time2));

% Low-pass filter SST
sstLP = NaN(size(sstSw1));
for i=1:length(lon)
    sstLP(i,:,:) = (pl66tn(squeeze(sstSw1(i,:,:)),1,240))';
end
%%
% Make a climatological annual cycle for each lat, lon pair
dn = datenum(time1); % convert to datenum
dv = datevec(time1); % Convert to datevec
yd = dn - datenum(dv(:,1),1,1) + 1; % Convert to yearday
% Vector to use for matching times to 6-hourly data
yhr = 1:1/4:(367-1/4);
% Vector to store climatology
sstSw0 = NaN(size(sstSw1));

% loop by values
for i=1:366*4
   times = ismembertol(yd,yhr(i)); % Logical of times in the 42-yr record that match the yearday and time
   mu = mean(sstSw1(:,:,times),3,'omitnan'); % take the mean along time dimension
   k = find(times); % Indices of the nonzero points in times
   for j = 1:length(k) % Loop through these times
       sstSw0(:,:,k(j)) = mu; % Assign 2D mean array to the time slice
   end
end 

% Low-pass filter climatological annual cycle
for l=1:length(lon) % loop along one spatial dimension since pl66 can only filter 2D arrays,
    % we can use either as long as time is the longer dimension
    sstSw0(l,:,:) = (pl66tn(squeeze(sstSw0(l,:,:)),1,240))'; % apply filter and assign to 2D longitude slice
end

%%
% Take difference between sstSw1 and sstSw0 to find SST'
sstSwA = sstSw1-sstSw0;

% Bandpass SST' by applying 6-month high-pass filter

% 6-month (half of a year) low-pass filter
hrs = hours(years(0.5)); % define the cutoff frequency in hours
lp6mo = nan(size(sstSwA)); % Initialize empty matrix to hold low-pass part
for m=1:length(lon) % loop along one spatial dimension since pl66 can only filter 2D arrays,
    % we can use either as long as time is the longer dimension
    lp6mo(m,:,:) = pl66tn(squeeze(sstSwA(m,:,:)),1,hrs)'; % Evaluate the low-pass filtered signal  and assign to 2D longitude slice
end
% Take high-pass part of signal
sstSwA = sstSwA-lp6mo;
% Replace one window-length with NaNs on each end
sstSwA(:,:,1:2*round(hrs))=NaN;
sstSwA(:,:,end-2*round(hrs):end)=NaN;

%%
% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(lat) max(lat)];
lonlim = [min(lon) max(lon)];
load coastlines

% Plot contours of SST' on world map for each date in summerDates
for n = 1:length(summerDates)
    t0 = ismembertol(time1,summerDates(n),hours(3)); % Finds which 6-hourly time is nearest to summerDates(n)
    % Since our tolerance is 3 hours the window of tolerance around each
    % summerDate is 6 hours wide, there is no possibility of having 2 t0's
    figure(n)
    worldmap(latlim,lonlim) % Map over Chile-Peru System
    plotm(coastlat,coastlon) % Adds coastlines
    [C,h] = contourm(Lat,Lon,sstSwA(:,:,t0)); % Contour of SST'
end


