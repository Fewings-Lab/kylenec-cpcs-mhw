% s02_1_anomaly_2016example_KMCsample.m 

% A script to load SST' and wind stress anomaly data for the Jan 2016
% example in Cooley et al. (2022) and then map SST' on 16 Jan 2016 and the
% mean of the wind stress anomaly vectors over 5-16 Jan 2016 (Figure 2).

% Version 20 Jan 2023 Kylene M Cooley
% Revised to share as sample of work. Removed intermediate plotting
% sections used while developing this script, and added context for data
% importing methods.

% Uses additional libraries and functions:
% landmaskcube - Kylene Cooley
% avgmap - Kylene Cooley
% read_era5_new - Larry O'Neill
% quivermc - Chad Greene
% Mapping Toolbox - MATLAB

%% Load anomaly data

% Define region to be displayed in maps
latlim = [-50 -15];                                                 % Southern Chilean to southern Peruvian coastline
lonlim = 360 + [-90 -70];                                           % Southeast Pacific to coast of South America

% Retrieve a subset of ERA5 unfiltered SST anomaly
SST = read_era5_new('sfc','daily','01/16/2016','01/16/2016',...     % SST anomaly
    {'sst'},lonlim,latlim,'anomaly',1);

% NOTE: Function read_era5_new() was used to import a subset of ERA5 data
% that was hosted on a remote OSU server. The ERA5 dataset there was
% updated regularly from the Climate Data Store. Data analysis was
% performed on the same OSU server.
% 'sfc' specifies the ocean surface level since ERA5 also provides some
% variables at different depths or elevations. 
% StartDate and EndDate are the same so that we only recieve one daily value. 
% Last argument is boolean 1 or 0 for true or false for anomaly or raw
% data. This and 'anomaly' is unnecessary when we don't want the anomaly.

% Define dates of ERA5 surface wind stress anomaly to retrieve
stdate = '01/05/2016';                                              % start date of example event
endate = '01/16/2016';                                              % end date of example event

% Retrieve a subset of unfiltered eastward (u) and northward (v) wind
% stress anomalies and the magnitude anomalies for surface wind stress
% vectors
uwindstr = read_era5_new('sfc','daily',stdate,endate,...            % zonal surface wind stress anomaly component
    {'ustr_sfc'},lonlim,latlim,'anomaly',1);   
vwindstr = read_era5_new('sfc','daily',stdate,endate,...            % meridional surface wind stress anomaly component
    {'vstr_sfc'},lonlim,latlim,'anomaly',1);   
windstrm = read_era5_new('sfc','daily',stdate,endate,...            % surface wind stress magnitude anomaly
    {'strm_sfc'},lonlim,latlim,'anomaly',1);   

%% Average daily surface wind stress anomaly vector components and magnitude for cmap
% The following means are evaluated along the time dimension so that the
% output is a 2D array of the corresponding quantity averaged at each
% lat-lon location over the twelve days 5 Jan to 16 Jan 2016.

ustr = mean(uwindstr.anomaly_ustr_sfc,3);                           % zonal component surface wind stress magnitude
vstr = mean(vwindstr.anomaly_vstr_sfc,3);                           % meridional component surface wind stress magnitude
strm = mean(windstrm.anomaly_strm_sfc,3);                           % surface wind stress vector magnitude

%% Load observed winds for climatology and confidence interval
% Retrieves the observed surface wind stress magnitude 01/05-01/16 for each
% year available from ERA5. The confidence interval on the climatology (or
% mean wind stress magnitude) for this period was calculated to compare
% with the value of the mean anomaly wind stress magnitude for
% 01/05/2016-01/16/2016. Mean wind stress magnitude anomalies (evaluated
% above) that are not outside the confidence interval for the climatology
% centered on zero are considered not significantly different from the
% climatology and replaced with NaN.

% Apply land mask to the wind stress magnitude 
strm = landmaskcube(strm,SST.lat,SST.lon);

% Initial parameters:
yrs = string(1980:1:2019);                                          % years to loop through below
% NOTE: These are the years that were included in our time series used to
% identifying events after the bandpass filter was applied.
stdy = '01/05';                                                     % start date as a string in format "mm/dd"
endy = '01/16';                                                     % end date as a string in format "mm/dd"
latlim = [-50 -15];                                                 % latitude limits of the subset region to retrieve (same as above)
lonlim = [270 290];                                                 % longitude limits of the subset region to retrieve

% Empty array to store wind stress magnitude values:
windavg = nan((range(latlim))*4+1,(range(lonlim))*4+1,length(yrs));
% NOTE: The size of this array is defined with previous knowledge that we
% are using a quarter-degree grid, so the total latitudinal and
% longitudinal range needs to be multiplied by 4 and 1 is added to each in
% order to have enough room to include the endpoints of the range.

% Loop to load the daily mean wind stress magnitude for 5-16 Jan of each
% year and evaluate the mean along the time dimension
for i = 1:length(yrs)
    % Start date and end date strings
    stdt = join([stdy,yrs(i)],"/");                                 % Joins stdy with year string to form a full date string with format "mm/dd/yyy"
    endt = join([endy,yrs(i)],"/");                                 % Joins endy with year string with format "mm/dd/yyyy"
    
    % Load 5-16 Jan from ERA5 for yrs(i)
    winds = read_era5_new('sfc','daily',stdt,endt,{'strm_sfc'},...
        lonlim,latlim);
    
    % Assign time average to a slice of windavg array for the corresponding year
    windavg(:,:,i) = squeeze(mean(winds.strm_sfc,3));
end

% Evaluate the confidence interval on the climatology (mean of windavg over time)
alpha = 0.05;                                                       % defines confidence level (1-alpha)*100%
sigma = std(windavg,1,3);                                           % sample estimate of the variance at each lat-lon location over all years
N = length(windavg(1,1,:));                                         % degrees of freedom
p = 1-(alpha/2);                                                    % x-value along distribution where upper tail is evaluated from
q_t = tinv(p,N-1);                                                  % upper tail of students-t distribution with N-1 effective degrees of freedom

delta_mu = (q_t/sqrt(N)).*sigma;                                    % uncertainty in climatology/half of the confidence interval

% Remove values of the mean wind stress magnitude anomaly that are in the
% interval [-delta_mu, +delta_mu] about zero (inclusive)
boo = delta_mu>=abs(strm);                                          % boolean where the above condition is true
strm_sig = strm;                                                    % assign same wind stress magnitude anomaly values to new statistically significant array
strm_sig(boo) = nan;                                                % ignore values where boo is true

% NOTE: Most of the wind stress magnitudes are very small compared to
% global max and min values in strm_sig. The number of points where
% strm_sig is larger than the value used as the cap is small (less than 10
% points).

% Cap wind stress magnitude values at +-0.15 Pa for plotting 
strm_sig(strm_sig>=0.15) = 0.15;                                    % upper value cap
strm_sig(strm_sig<= -0.15) = -0.15;                                 % lower value floor

%% Create subplots with SST anomaly and then wind stress anomaly and titles are sub-figure labels
% Maps the SST anomaly and statistically significant surface wind stress
% anomalies for the Jan 2016 example figure (Figure 2) as subplots in the
% same figure.

% Define lat and lon arrays for the wind stress anomaly vectors 
Lat = SST.lat*ones(length(SST.lon),1)';                             % 2D lat array
Lat = Lat(1:15:end,1:15:end);                                       % subset of lat array for wind stress anomaly vectors
Lon = ones(length(SST.lat),1)*SST.lon';                             % 2D lon array
Lon = Lon(1:15:end,1:15:end);                                       % subset of lon array for wind stress anomaly vectors

% Load MATLAB coastline .mat file. Mapping toolbox must be installed in packages. 
load coastlines.mat

% Apply land-sea-mask to wind stress anomaly vector arrays 
ustr = landmaskcube(ustr,SST.lat,SST.lon);                          % apply land mask to zonal wind stress anomaly component
vstr = landmaskcube(vstr,SST.lat,SST.lon);                          % apply land mask to meridional wind stress anomaly component
U = ustr(1:15:end,1:15:end);                                        % zonal subset for wind stress anomaly vectors
V = vstr(1:15:end,1:15:end);                                        % meridional subset for wind stress anomaly vectors

% Open figure
figure()

% Figure 2a
ax1 = subplot(1,2,1);                                               % create subplot axes
[Csst,hsst,csst] = avgmap(SST.anomaly_sst,SST.lat,SST.lon,...       % plot colormap of SST anomaly in ax1
    'a) 16 Jan 2016',"$\textsf{SST anomaly [}^\circ\textsf{C]}$",...
    'balance',ax1); 
title1 = ax1.Title;                                                 % assign title text object to title1
title1.FontSize = 8;                                                % set title text font size
title1.Units = 'Normalize';                                         % use normalized figure units for title location (horizontal and vertical lengths of axes are on interval [0, 1])
title1.Position(1) = -0.2;                                          % shift lower left corner of title location to left
ax1.TitleFontWeight = 'normal';                                     % set title font weight
ax1.TitleHorizontalAlignment = 'left';                              % set title alignment (in text box)
setm(ax1,'FontSize',8)                                              % set font size for everything in axes ax1
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');                       % plot white dot over Punta Lavapie, Chile
PL1.Children.ZData = 4;                                             % set z-level of dot so that it lies over coastline and data
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');  % plot black dot in center of white dot
PL2.Children.ZData = 5;                                             % set z-level of dot to sit on top of white dot
t3 = text(610000,-4250000,'PL','Color',[0 0 0],...                  % add label for Punta Lavapie dot
    'HorizontalAlignment','left','FontSize',8); 
t1 = text(-332137,-2840660,'warmer','Color',[1 1 1],...             % add text over warmer SST' region in map
    'HorizontalAlignment','center','FontSize',8); 

% Figure 2b
ax2 = subplot(1,2,2);                                               % go to second subplot axes
[Cstrm,hstrm,cstrm] = avgmap(strm_sig,SST.lat,SST.lon,...           % plot colormap of wind stress magnitude anomaly in ax2
    'b) 5-16 Jan 2016', "$\textsf{Wind stress anomaly [Pa]}$",...
    'balance',ax2);           
cstrm.Ticks = [cstrm.Ticks 0.15];                                   % append max wind stress anomaly value to top of colorbar
quivermc(Lat,Lon,U,V,'units','N m^{-2}')                            % overlay wind stress anomaly vectors
title2 = ax2.Title;                                                 % assign title text object to tit2
title2.FontSize = 8;                                                % set title font size
title2.Units = 'Normalize';                                         % convert title location units to normalized units w.r.t. ax2
title2.Position(1) = -0.2;                                          % shift lower left corner of title location to left
ax2.TitleFontWeight = 'normal';                                     % set title font weight (removes bold)
ax2.TitleHorizontalAlignment = 'left';                              % set title alignment 
setm(ax2,'FontSize',8)                                              % set font size for all of ax2
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');                       % plot white dot over Punta Lavapie
PL1.Children.ZData = 4;                                             % set white dot zdata to overlie everything else
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');  % plot black dot over white dot
PL2.Children.ZData = 5;                                             % set black dot zdata to sit on top of white dot
t4 = text(610000,-4250000,'PL','Color',[0 0 0],...                  % add label for Punta Lavapie dot
    'HorizontalAlignment','left','FontSize',8); 
t2 = text(-332137,-2840660,'weaker','Color',[1 1 1],...             % add text over weaker wind stress magnitude anomaly region
    'HorizontalAlignment','center','FontSize',8); 

% Set figure position on computer screen and print position to full
% text-width and desired height
set(gcf,'Units','centimeters','Position',[0 0 19 11],...            % sets position to 19 cm wide by 11 cm tall in lower left corner
    'PaperUnits','centimeters','PaperPosition',[0 0 19 11]) 
% NOTE: Position vector is defined as [LowerLeftCornerX, LowerLeftCornerY, Width, Height]

% Set colorbar text size to 8 pt 
cstrm.FontSize = 8;                                                 % wind stress magnitude colorbar
csst.FontSize = 8;                                                  % SST colorbar

% Export figure to PDF with 2020 pixel resolution (required by AGU for
% paper size defined above)
% exportgraphics(gcf,'2016event_anomaly_v10.pdf','Resolution',2020)
% exportgraphics(gcf,'2016event_anomaly_v10.png','Resolution',2020)