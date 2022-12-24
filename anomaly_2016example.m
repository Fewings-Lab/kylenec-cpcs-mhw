% anomaly_2016example.m 

% A script to load SST' and wind stress anomaly data for the Jan 2016
% example and then map SST' on 16 Jan 2016 and the mean of the wind stress
% anomaly vectors over 5-16 Jan 2016. The section plotting these as two
% subplots in a single figure yields Figure 2 in Cooley et al. (2022)

% Version 8/26/2022 Kylene M Cooley

%% Load anomaly data

% Define region to be displayed in maps
latlim = [-50 -15];         % Southern Chilean to southern Peruvian coastline
lonlim = 360 + [-90 -70];   % Southeast Pacific to coast of South America

% Retrieve a subset of ERA5 unfiltered SST anomaly
SST = read_era5_new('sfc','daily','01/16/2016','01/16/2016',{'sst'},lonlim,latlim,'anomaly',1); % SST anomaly
% NOTE: 'sfc' specifies the ocean surface level since ERA5 also provides some
% variables at different depths or elevations. StartDate and EndDate are
% the same so that we only recieve one daily value. Last argument is
% boolean 1 or 0 for true or false for anomaly or measured value. This and
% 'anomaly' is unnecessary when we don't want the anomaly.

% Define dates of ERA5 unfiltered surface wind stress anomaly to retrieve
stdate = '01/05/2016';
endate = '01/16/2016';
% NOTE: Figure 2b in Cooley et al. (2022) is an average of wind stress over
% this period in comparsion to the average of climatological values during
% these dates in January over the years available from ERA5. [This phrasing
% doesn't make a ton of sense, need to revise]

% Retrieve a subset of unfiltered eastward (u) and northward (v) wind
% stress anomalies and the magnitude anomalies for surface wind stress
% vectors
uwindstr = read_era5_new('sfc','daily',stdate,endate,{'ustr_sfc'},lonlim,latlim,'anomaly',1);   % zonal surface wind stress anomaly component
vwindstr = read_era5_new('sfc','daily',stdate,endate,{'vstr_sfc'},lonlim,latlim,'anomaly',1);   % meridional surface wind stress anomaly component
windstrm = read_era5_new('sfc','daily',stdate,endate,{'strm_sfc'},lonlim,latlim,'anomaly',1);   % surface wind stress magnitude anomaly

%% Map of SST' (Figure 2a)
% The figure plotted in this section has been included in the last section
% as the first subplot in that figure. avgmap.m now has an optional argument to
% call an axes object to plot the data on. 

% avgmap(SST.anomaly_sst,SST.lat,SST.lon,' ',"$SST'$ [$^\circ$C]",'balance')

%% Average daily surface wind stress anomaly vector components and magnitude for cmap
% The following means are evaluated along the time dimension so that the
% output is a 2D array of the corresponding quantity averaged at each
% lat-lon location over the twelve days 5 Jan to 16 Jan 2016.

ustr = mean(uwindstr.anomaly_ustr_sfc,3); % zonal component surface wind stress magnitude
vstr = mean(vwindstr.anomaly_vstr_sfc,3); % meridional component surface wind stress magnitude
strm = mean(windstrm.anomaly_strm_sfc,3); % surface wind stress vector magnitude

%% Load observed winds for climatology and confidence interval
% This section retrieves the observed surface wind stress magnitude
% 01/05-01/16 for each year available from ERA5. From this we evaluated the
% confidence interval on the climatology (or mean wind stress magnitude
% over time) to compare with the value of the mean anomaly wind stress
% magnitude for 01/05/2016-01/16/2022. Mean wind stress magnitude anomalies
% (evaluated above) that are not outside the confidence interval for the
% climatology centered on zero are considered not significantly different
% from the climatological surface wind stress magnitude at each lat-lon
% pair location.

% Apply land mask to the wind stress magnitude 
strm = landmaskcube(strm,SST.lat,SST.lon);

% Initial parameters:
yrs = string(1980:1:2019); % years to loop through below
% NOTE: We defined which years to use in this study based on the years that were included
% in our time series for finding events after the bandpass filter was
% applied. Another array containing strings of years in each cell could be
% used instead.
stdy = '01/05';             % start date as a string in format "mm/dd"
endy = '01/16';             % end date as a string in format "mm/dd"
latlim = [-50 -15];         % latitude limits of the subset region to retrieve (same as above)
lonlim = [270 290];         % longitude limits of the subset region to retrieve

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
    stdt = join([stdy,yrs(i)],"/");             % Joins stdy with year string to form a full date string with format "mm/dd/yyy"
    endt = join([endy,yrs(i)],"/");             % Joins endy with year string with format "mm/dd/yyyy"
    
    % Load 5-16 Jan from ERA5 for yrs(i)
    winds = read_era5_new('sfc','daily',stdt,endt,{'strm_sfc'},lonlim,latlim);
    
    % Assign time average to a slice of windavg array for the corresponding year
    windavg(:,:,i) = squeeze(mean(winds.strm_sfc,3));
end

% Evaluate the confidence interval on the climatology (mean of windavg over time)
alpha = 0.05;                       % defines confidence level (1-alpha)*100%
sigma = std(windavg,1,3);           % sample estimate of the variance at each lat-lon location over all years
N = length(windavg(1,1,:));         % degrees of freedom
p = 1-(alpha/2);                    % x-value along distribution where upper tail is evaluated from
q_t = tinv(p,N-1);                  % upper tail of students-t distribution with N-1 effective degrees of freedom

delta_mu = (q_t/sqrt(N)).*sigma;    % uncertainty in climatology/half of the confidence interval

% Remove values of the mean wind stress magnitude anomaly that are in the
% interval [-delta_mu, +delta_mu] about zero (inclusive)
boo = delta_mu>=abs(strm);          % boolean where the above condition is true
strm_sig = strm;                    % assign same wind stress magnitude anomaly values to new statistically significant array
strm_sig(boo) = nan;                % ignore values where boo is true

% NOTE: most wind stress magnitudes are very small compared to global max and min
% values in strm_sig, while the number of points where strm_sig is
% uncharacteristically large is small (less than 10)

% Cap wind stress magnitude values at +-0.15 Pa for plotting 
strm_sig(strm_sig>=0.15) = 0.15;    % upper value cap
strm_sig(strm_sig<= -0.15) = -0.15; % lower value floor

%% Map average surface wind stress anomaly arrows and magnitude colormap (Figure 2b)
% The figure plotted in this section has been included in the last section
% as the second subplot in that figure. 

% % Map wind stress magnitude with colormap
% [Cstrm,hstrm,cstrm] = avgmap(strm_sig,SST.lat,SST.lon,' ',"$|\vec{\tau}|'$ [Pa]",'balance'); % [Csst,hsst,csst]=[map axes, patch object, colorbar object]
% cstrm.Ticks = [cstrm.Ticks 0.15];                                          % appends max wind stress magnitude to top of colorbar
% 
% % Create lat and lon arrays for plotting wind stress vectors
% Lat = SST.lat*ones(length(SST.lon),1)';                                    % latitude array with same size as wind stress arrays
% Lat = Lat(1:15:end,1:15:end);                                              % subset of latitudes for wind vectors
% Lon = ones(length(SST.lat),1)*SST.lon';                                    % longitude array with same size as wind stress arrays
% Lon = Lon(1:15:end,1:15:end);                                              % subset of longitudes for wind vectors
% 
% % Apply land-sea-mask to meridional and zonal wind stress anomaly vector components 
% ustr = landmaskcube(ustr,SST.lat,SST.lon);                                 % zonal wind stress anomaly component
% vstr = landmaskcube(vstr,SST.lat,SST.lon);                                 % meridional wind stress anomaly component
% U = ustr(1:15:end,1:15:end);                                               % subset of u for vectors to be mapped
% V = vstr(1:15:end,1:15:end);                                               % subset of v for vectors to be mapped
% 
% % Add wind stress anomaly vectors to wind stress anomaly magnitude color map
% h = quivermc(Lat,Lon,U,V,'units','N m^{-2}'); 

%% Create subplots with SST anomaly and then wind stress anomaly and titles are sub-figure names (Figure 2)
% Maps the SST anomaly and statistically significant surface wind stress
% anomalies for the Jan 2016 example figure (Figure 2) as subplots in the
% same figure.

% Define lat and lon arrays for the wind stress anomaly vectors 
Lat = SST.lat*ones(length(SST.lon),1)';                         % 2D lat array
Lat = Lat(1:15:end,1:15:end);                                   % subset of lat array for wind stress anomaly vectors
Lon = ones(length(SST.lat),1)*SST.lon';                         % 2D lon array
Lon = Lon(1:15:end,1:15:end);                                   % subset of lon array for wind stress anomaly vectors

% Load MATLAB coastline .mat file. Mapping toolbox must be installed in packages 
load coastlines.mat

% Apply land-sea-mask to wind stress anomaly vector arrays 
ustr = landmaskcube(ustr,SST.lat,SST.lon);                      % apply land mask to zonal wind stress anomaly component
vstr = landmaskcube(vstr,SST.lat,SST.lon);                      % apply land mask to meridional wind stress anomaly component
U = ustr(1:15:end,1:15:end);                                    % zonal subset for wind stress anomaly vectors
V = vstr(1:15:end,1:15:end);                                    % meridional subset for wind stress anomaly vectors

% Open figure
figure()

% Figure 2a
ax1 = subplot(1,2,1);                                               % create subplot axes
[Csst,hsst,csst] = avgmap(SST.anomaly_sst,SST.lat,SST.lon,'a) 16 Jan 2016',"$\textsf{SST anomaly [}^\circ\textsf{C]}$",'balance',ax1); % plot colormap of SST anomaly in ax1
tit1 = ax1.Title;                                                   % assign title text object to tit1
tit1.FontSize = 8;                                                  % set title text font size
tit1.Units = 'Normalize';                                           % use normalized figure units for title location (horizontal and vertical lengths of axes are on interval [0, 1])
tit1.Position(1) = -0.2;                                            % shift lower left corner of title location to left
ax1.TitleFontWeight = 'normal';                                     % set title font weight
ax1.TitleHorizontalAlignment = 'left';                              % set title alignment (in text box)
setm(ax1,'FontSize',8)                                              % set font size for everything in axes ax1
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');                       % plot white dot over Punta Lavapie, Chile
PL1.Children.ZData = 4;                                             % set z-level of dot so that it lies over coastline and data
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');  % plot black dot in center of white dot
PL2.Children.ZData = 5;                                             % set z-level of dot to sit on top of white dot
t3 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8); % add label for Punta Lavapie dot
t1 = text(-332137,-2840660,'warmer','Color',[1 1 1],'HorizontalAlignment','center','FontSize',8); % add text over warmer SST' region in map

% Figure 2b
ax2 = subplot(1,2,2);                                               % go to second subplot axes
[Cstrm,hstrm,cstrm] = avgmap(strm_sig,SST.lat,SST.lon,'b) 5-16 Jan 2016',"$\textsf{Wind stress anomaly [Pa]}$",'balance',ax2); % plot colormap of wind stress magnitude anomaly in ax2
cstrm.Ticks = [cstrm.Ticks 0.15];                                   % append max wind stress anomaly value to top of colorbar
quivermc(Lat,Lon,U,V,'units','N m^{-2}')                            % overlay wind stress anomaly vectors
%fillm(coastlat,coastlon,[0.7 0.7 0.7])                             % Adds gray landmass - included in avgmap function
tit2 = ax2.Title;                                                   % assign title text object to tit2
tit2.FontSize = 8;                                                  % set title font size
tit2.Units = 'Normalize';                                           % convert title location units to normalized units w.r.t. ax2
tit2.Position(1) = -0.2;                                            % shift lower left corner of title location to left
ax2.TitleFontWeight = 'normal';                                     % set title font weight (removes bold)
ax2.TitleHorizontalAlignment = 'left';                              % set title alignment 
setm(ax2,'FontSize',8)                                              % set font size for all of ax2
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');                       % plot white dot over Punta Lavapie
PL1.Children.ZData = 4;                                             % set white dot zdata to overlie everything else
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');  % plot black dot over white dot
PL2.Children.ZData = 5;                                             % set black dot zdata to sit on top of white dot
t4 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8); % add label for Punta Lavapie dot
t2 = text(-332137,-2840660,'weaker','Color',[1 1 1],'HorizontalAlignment','center','FontSize',8); % add text over weaker wind stress magnitude anomaly region

% Set figure position on computer screen and print position to full
% text-width and desired height
set(gcf,'Units','centimeters','Position',[0 0 19 11],'PaperUnits','centimeters','PaperPosition',[0 0 19 11]) % sets position to 19 cm wide by 11 cm tall in lower left corner
% NOTE: Position vector is defined as [LowerLeftCornerX, LowerLeftCornerY, Width, Height]

% Set colorbar text size to 8 pt 
cstrm.FontSize = 8;                     % wind stress magnitude colorbar
csst.FontSize = 8;                      % SST colorbar

% Export figure to PDF with 2020 pixel resolution (required by AGU for
% paper size defined above)
% exportgraphics(gcf,'2016event_anomaly_v10.pdf','Resolution',2020)
% exportgraphics(gcf,'2016event_anomaly_v10.png','Resolution',2020)