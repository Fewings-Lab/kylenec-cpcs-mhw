% SST_anomaly_maps.m
% Kylene Cooley 
% 29 Jan 2021 
% A script to maps the spatial extent of strong warm SST' at times we found during the
% strong upwelling season (Dec-Feb)

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
% time2 = interp1(1:6:6*length(time1),time1,1:6*length(time1)); % interpolates 6-hrly time to hourly, but not exactly on the hour
% time2 = datetime(time2,'ConvertFrom','datenum'); % convert to datetimes
% time2 = dateshift(time2, 'start', 'hour', 'nearest'); % shifts datetimes to nearest hour
time1 = datetime(time1,'ConvertFrom','datenum'); % convert non-interpolated time to datetimes

% % Save the relevant variables to a .mat file
% save('ChileCoast_SST_6H.mat','lat','lon','sstSw1','time1','summerDates')

% [Ignore this because it requires too much RAM:]
% Linearly interpolate 6-hrly SST (constant daily values) onto hourly grid
% sstSw2 = interp3(lat,lon,time1,sstSw1,lat,lon,datenum(time2));


% Import variables from ChileCost_SST_6H.mat
% load('ChileCoast_SST_6H.mat')
% Contains: sstSw1 in degC, lat [-45 -20], lon [-90 -70], time1 which is a
% 6-hourly datetime vector, and summerDates which is a datetime vector of
% times we want to plot at the end of this script

% Low-pass filter SST with a 10-day pl66 filter
sstLP = NaN(size(sstSw1)); % This will be the 10 low-pass filtered SST swath
for i=1:length(lon) % Looping through longitude slices because the pl66 filter only handles 2D matrices
    sstLP(i,:,:) = (pl66tn(squeeze(sstSw1(i,:,:)),6,240))'; % Don't forget that for 6-hourly data, dt=6
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

% loop by yearday w/ time
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
    sstSw0(l,:,:) = (pl66tn(squeeze(sstSw0(l,:,:)),6,240))'; % apply filter and assign to 2D longitude slice
end

%%
% Take difference between sstLP and sstSw0 to find SST'
sstSwA = sstLP-sstSw0;

% Bandpass SST' by applying 6-month high-pass filter

% 6-month (half of a year) low-pass filter
hrs = hours(years(0.5)); % define the cutoff frequency in hours
lp6mo = nan(size(sstSwA)); % Initialize empty matrix to hold low-pass part
for m=1:length(lon) % loop along one spatial dimension since pl66 can only filter 2D arrays,
    % we can use either as long as time is the longer dimension
    lp6mo(m,:,:) = pl66tn(squeeze(sstSwA(m,:,:)),6,hrs)'; % Evaluate the low-pass filtered signal  and assign to 2D longitude slice
end
% Take high-pass part of signal
sstSwA = sstSwA-lp6mo;
% Replace one window-length with NaNs on each end
sstSwA(:,:,1:2*round(hrs/6))=NaN;
sstSwA(:,:,end-2*round(hrs/6):end)=NaN;

%%
% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(lat) max(lat)];
lonlim = [min(lon) max(lon)]; % works for negative longitude, but will have to switch min/max if lon is in [0 360]
load coastlines
Lat = ones(length(lon),1)*lat';
Lon = lon*ones(1,length(lat));
clim = [min(sstSwA,[],'all') max(sstSwA,[],'all')];

A = datenum(time1);
% Plot contours of SST' on world map for each date in summerDates
for n = 1:length(summerDates)
    B = datenum(summerDates(n));
    tol = datenum(hours(3))/max(abs([A(:);B(:)]));
    
    t0 = ismembertol(A,B,tol); % Finds which 6-hourly time is nearest to summerDates(n)
    % [Here is where I went wrong and this is what I thought would happen:] 
    % Since our tolerance is 3 hours, the window of tolerance around each
    % summerDate is 6 hours wide, there is no possibility of having 2 t0's
    % [What actually happened: the logical t0 was true for every time and
    % that should not have happened with the tolerance so small compared to
    % the dates]
    if n<=16
        figure(1)
        subplot(4,4,n)
        worldmap(latlim,lonlim) % Map over Chile-Peru System
        plotm(coastlat,coastlon) % Adds coastlines
        if sum(t0)>1
            [C,h] = contourm(Lat,Lon,mean(sstSwA(:,:,t0),3,'omitnan'),'Fill','on'); % Contour of SST'
            caxis(clim)
        else
            [C,h] = contourm(Lat,Lon,sstSwA(:,:,t0),'Fill','on'); % Contour of SST'
            caxis(clim)
        end
        title(datestr(summerDates(n)))
    elseif n>=17
        figure(2)
        subplot(4,4,n-16)
        worldmap(latlim,lonlim) % Map over Chile-Peru System
        plotm(coastlat,coastlon) % Adds coastlines
        if sum(t0)>1
            [C,h] = contourm(Lat,Lon,mean(sstSwA(:,:,t0),3,'omitnan'),'Fill','on'); % Contour of SST'
            caxis(clim)
        else
            [C,h] = contourm(Lat,Lon,sstSwA(:,:,t0),'Fill','on'); % Contour of SST'
            caxis(clim)
        end
        title(datestr(summerDates(n)))
    end
end

figure(1)
sgtitle("SST'")
hp4 = get(subplot(4,4,16),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1])
figure(2)
sgtitle("SST'")
hp4 = get(subplot(4,4,16),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*4.1])


