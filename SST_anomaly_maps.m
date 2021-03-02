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
% Vector to use for matching times 
yhr = 1:1/4:(367-1/4);
% Vector to store climatology
sstSw0 = NaN(size(sstSw1));
for j = 1:length(lon)
    for k = 1:length(lat)
        foo2 = NaN(size(yd));
        % loop by values
        for i=1:366*4
           values = ismembertol(yd,yhr(i));
           mu = nanmean(sstSw1(j,k,values));
           foo2(values) = mu;
        end 
       
        
    end
end

% Filter the working variable
foo4 = pl66tn(foo3,1,240); % apply 10-day filter

%%
% Coastline data set and coordinate limits around N. Pacific Basin:
latlim = [10 70];
lonlim = [120 260];
load coastlines

% Plot contours of SST' on world map for each date in summerDates
figure(1)
worldmap(latlim,lonlim) % Map around N. Pacific
plotm(coastlat,coastlon) % Adds coastlines
[C,h] = contourm(Lat,Lon,sstA); % Contour of SLP levels for every 2 mbar