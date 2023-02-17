% s01_0_summer_SST_wind_climatology_CPCS.m

% A script to map summer (DJF) average SST' and wind stress as a for the
% Chile-Peru Current System between 15 S and 50 S. The section plotting
% these on the same map axes yields Figure 1 in Cooley et al. (2022).

% Version 9/2/2022 Kylene M Cooley

% NOTE: Consider changing name to something more descriptive of resulting
% figure

%% Define map region and time period to average summer months over, then collect SST data

% Initial variables defining our research area and time period
% NOTE: MATLAB requires monotonically increasing endpoints in a range
latlim = [-50 -15];             % latitude limits of 50 S to 15 S
lonlim = 360 + [-90 -70];       % longitude limits of 90 W to 70 W
stdate = '01/01/1979';          % start date is first date available from ERA5
endate = '12/31/2020';          % end date of data to retrieve (inclusive)
% NOTE: Dates were chosen to include as much data available in summer 2021
% as possible while retrieving only complete years.

% Collect non-anomaly SST in the study region and anomaly
go = isfile('SST_era5_v1.mat'); % boolean to determine whether .mat file with the data already exists in current working directory

if go == 0                      % if file containing data does not exist in current working directory...
    % Collect unfiltered eastward (u) and northward (v) wind stress fields
    % EXCEPT THERE'S NO COMMAND TO COLLECT WIND HERE, WHAT HAPPENED???
    % Collect unfiltered SST and SST anomaly in the area defined by latlim and lonlim above
    % between stdate and endate.
    sst = read_era5_new('sfc','daily',stdate,endate,{'sst'},lonlim,latlim);                         % daily SST, highest observation frequency available from ERA5
    sst_anomaly = read_era5_new('sfc','daily',stdate,endate,{'sst'},lonlim,latlim,'anomaly',1);     % daily SST anomaly
    save('SST_era5_v1.mat','lonlim','latlim','sst','sst_anomaly','-v7.3')                           % saves variables in workspace listed to compressed (version 7.3) .mat file
end

%% Evaluate the summer mean of SST

% Load the existing .mat file with SST and SST anomaly data to workspace
load SST_era5_v1.mat 

% Extract double type arrays from struct object with dot indexing
SST = sst.sst-273.15;                           % double type array with SST converted to degC from Kelvin
lat = sst.lat;                                  % double type vector of latitude values, dim 1 of SST
lon = sst.lon;                                  % double type vector of longitude values, dim 2 of SST
time = sst.time;                                % double type vector of observation times, dim 3 of SST

% Remove values not observed during summer (DJF)
[~,M,~,~,~,~] = datevec(time);                  % convert time datenums to date vector with 6 columns, and assign month column to M
keep = M > 11 | M < 3;                          % logical array same length as M, true (=1) when month number is 12, 1, or 2
summer_sst = SST(:,:,keep);                     % use logical keep along time dim to assign only summer SST values to new variable

% Apply land-sea-mask in case data contained some values over land
summer_sst = landmaskcube(summer_sst,lat,lon);  % function landmaskcube(3dArray,latVector,lonVector) applies mask array by looping along dim 3

% Evaluate mean of summer SST along time dimension (dim 3)
avg_sst = mean(summer_sst,3);

%% Summer mean of northward and eastward surface wind stress vector components

% Load existing eastward (u) and northward (v) wind stress data structs to
% workspace
load windstress_era5_v2.mat

% Extract each double type array from the corresponding variable struct
vstr = vwindstr.vstr_sfc;                           % double array with northward wind surface stress
ustr = uwindstr.ustr_sfc;                           % double array with eastward wind surface stress

% Isolate summer wind stress values with same 'keep' logical vector along time dimension
summer_vstr = vstr(:,:,keep);                       % 3D summer northward surface wind stress array
summer_ustr = ustr(:,:,keep);                       % 3D summer eastward surface wind stress array

% Apply land-sea-mask to wind stress vector components as above
summer_vstr = landmaskcube(summer_vstr,lat,lon);    % summer northward surface wind stress
summer_ustr = landmaskcube(summer_ustr,lat,lon);    % summer eastward surface wind stress

% Evaluate mean of summer surface wind stress vector components along time
% dimension (dim 3)
avg_vstr = mean(summer_vstr,3);                     % 2D summer average northward surface wind stress, avg_vstr(lat,lon)
avg_ustr = mean(summer_ustr,3);                     % 2D summer average eastward surface wind stress, avg_ustr(lat,lon) 

%% Plot summer mean SST map
% Additional map elements in Figure 1 of Cooley et al. (2022) are added to
% this map in subsequent sections.

% Set different argument variables for mapping function
titlestr = "Average summer SST and surface wind stress";    % define the title
clabel = "$\textsf{SST [}^{\circ}\textsf{C]}$";             % define the color bar label
cmap = 'thermal';                                           % set color map to use from cmocean package

% Open figure window and plot average SST on map axes
[~,~,~] = avgmap(avg_sst,lat,lon,[],clabel,cmap);           % function avgmap(2dArray,lat,lon,title,colorbarLabel,colormap)

%% Overlay surface wind stress vectors

Lat = lat*ones(length(lon),1)';
Lat = Lat(1:15:end,1:15:end);
Lon = ones(length(lat),1)*lon';
Lon = Lon(1:15:end,1:15:end);

% Lat = [Lat nan(size(Lat,1),1)];
% Lon = [Lon nan(size(Lon,1),1)];
% Lat(end,end) = -47.5;
% Lon(end,end) = -72.5+360;

U = avg_ustr(1:15:end,1:15:end);
V = avg_vstr(1:15:end,1:15:end);

% U = [U nan(size(U,1),1)];
% V = [V nan(size(V,1),1)];
% U(end,end) = 0.05;
% V(end,end) = 0;

% quivers of wind speed are plotted over SST with quivermc b/c quiverm can
% only accept changes in degree lat or degree lon as u or v and switches
% the order of u and v from usual geophysical use. quivermc can use
% absolute u,v values but these are plotted over the lat,lon locations.
% Lat, lon, u, or v inputs below must be grids with 2 dimensions > 1.
h = quivermc(Lat,Lon,U,V,'units','N m^{-2}'); 

%% Overlay nearshore box 
ptLat = [-35.5 -29]; % location of time series points
ptLon = [-72.75 -77.5];
h1 = patchm([ptLat(1) ptLat(1)-1 ptLat(1)-1 ptLat(1)],[ptLon(1)-(3/5)  ptLon(1)-(3/5)  ptLon(1)-(8/5) ptLon(1)-(8/5)],'EdgeColor','g','FaceColor','none','LineWidth',1.5);

% set(gcf,'PaperPosition',[0 0 5 6])

% % Plotting reference quiver arrows
% testlon_a = [-72.5,-72.5;-86,-86]+360;
% testlat_a = [-47.5,-25+1;-25,-25+1];
% u1 = [0.05,0;0,0];u2 = [0.05,0;0,0];u3 = [0,0;0,0];
% v1 = [0.05,0;0,0];v2 = [0,0;0,0];v3 = [0.05,0;0,0];
% % arrows cross over the lat, lon location
% quivermc(testlat_a,testlon_a,u1,v1);
% quivermc(testlat_a,testlon_a,u2,v2) 
% quivermc(testlat_a,testlon_a,u3,v3)
% 
% testlon_b = [-80,-80,-80]+360;
% testlat_b = [-45,-45,-45];
% % test whether quiverm distorts arrows with latitude
% g = quivermc(testlat_b,testlon_b,[3,0,3],[0,3,3],'g');


% saveas(gcf,'avg_SST_windstr_summer_withBox_v7.png')

%% Overlay offshore box for appendix version 
% about no influence of upwelling filaments

h2 = patchm([ptLat(2)+0.5 ptLat(2)-0.5 ptLat(2)-0.5 ptLat(2)+0.5],[ptLon(2)+0.5  ptLon(2)+0.5  ptLon(2)-0.5 ptLon(2)-0.5],'EdgeColor','c','FaceColor','none','LineWidth',1.5);
set(gcf,'Units','centimeters','Position',[0 0 9.5 11],'PaperUnits','centimeters','PaperPosition',[0 0 9.5 11])
set(gca,'FontSize',8)
setm(gca,'FontSize',8)
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');
PL1.Children.ZData = 4;
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');
PL2.Children.ZData = 5;
t1 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8);

% exportgraphics(gcf,'avg_SST_windstr_summer_offshoreBox_v9.pdf','Resolution',1600)
%% 
% lat_test = [-45.*ones(3,1) -26.*ones(3,1) NaN(3,1)];
% lon_test = [(-80.*ones(3,2)+360) NaN(3,1)];
% U_test = [1 1 0; 1 1 0; 0 0 0];
% V_test = [0 0 0; 1 1 0; 1 1 0];
% 
% h = quivermc(lat_test,lon_test,U_test,V_test); 
%% TUTORIAL: Determining if lengths of quiver arrows on a map are scaled with latitude
% When plotting arrows on a map, it is preferable to have [TBH I'M NOT SURE
% IF THEY SHOULD OR SHOULDN'T SCALE THE LENGTH WITH LATITUDE. IF THEY DO,
% WE CAN'T COMPARE ARROWS AT DIFFERENT LATITUDES. I THINK IT IS MORE LIKE
% DON'T SCALE WITH LATITUDE BUT KEEP THE ANGLES CONSISTENT?]

% Plotting reference quiver arrows
% testlon_a = [-80,-80;-86,-86]+360;
% testlat_a = [-26,-45;-26,-45];
% u1 = [1,0;0,0];u2 = [1,0;0,0];u3 = [0,0;0,0];
% v1 = [1,0;0,0];v2 = [0,0;0,0];v3 = [1,0;0,0];
% % arrows cross over the lat, lon location
% quivermc(testlat_a,testlon_a,u1,v1);
% quivermc(testlat_a,testlon_a,u2,v2) 
% quivermc(testlat_a,testlon_a,u3,v3)