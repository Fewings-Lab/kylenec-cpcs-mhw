% SST_analysis.m 
% Kylene Cooley 
% 29 Jan 2021 
% A script to carry out analysis of SST anomalies from a netcdf-format
% multi-file dataset

%% Load dataset from netcdf format
finfo = ncinfo('Data/Constitucion-SST-2018-2020.nc');
data1 = squeeze(ncread('Data/Constitucion-SST-1979-1991.nc','sst'));
data2 = squeeze(ncread('Data/Constitucion-SST-1992-2004.nc','sst'));
data3 = squeeze(ncread('Data/Constitucion-SST-2005-2017.nc','sst'));
data4 = ncread('Data/Constitucion-SST-2018-2020.nc','sst');
data4 = squeeze(data4(:,:,1,:));
lat = ncread('Data/Constitucion-SST-2018-2020.nc','latitude');
lon = ncread('Data/Constitucion-SST-2018-2020.nc','longitude');
time1 = ncread('Data/Constitucion-SST-1979-1991.nc','time');
time2 = ncread('Data/Constitucion-SST-1992-2004.nc','time');
time3 = ncread('Data/Constitucion-SST-2005-2017.nc','time');
time4 = ncread('Data/Constitucion-SST-2018-2020.nc','time');

% Combine like variables from different files
sst = cat(1,data1,data2,data3,data4);
time = cat(1,time1,time2,time3,time4);

% get datenum vector
tf = (length(time)/24)+1;
dn = datenum(1979,1,1:1/24:tf).';
dn = dn(1:end-1);

% Save variables as .mat file
save('ptLavapie_SST.mat','lat','lon','sst','time','dn')

%% Load data from .mat file
clear
load('ptLavapie_SST.mat')
sst = sst-273.15; % Convert sst to degC

%% Use pl66 to low-pass filter sst
% with dt=1 (same as default) and T=168 hrs (7 days)
sstf = pl66tn(sst,1,168);

%% Obtain climatology
% Convert datenums to date vector [yy mm dd hh mm ss]
dv = datevec(dn); 

% Calculate yearday by subtracting the datenum of Jan 1 of that year 
% (and add 1 so the first day of the year is yearday 1 not 0)
yd = dn - datenum(dv(:,1),1,1) + 1; 

% Vectors to use for matching times and storing climatology
yhr = 1:1/24:(367-1/24);
sst0 = NaN(length(yd),1);

% loop by values
for i=1:366*24
    values = ismembertol(yd,yhr(i));
    mu = nanmean(sstf(values));
    sst0(values) = mu;
end 
    
