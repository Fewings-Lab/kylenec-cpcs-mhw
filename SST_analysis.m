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
SST = timeseries(sst,datestr(dn)); % Make sst timeseries object
% % Plot raw data
% figure(1)
% plot(sst)
% title('ERA-5 SST at -35.51, -72.77')
% xlabel('Date')
% xtickformat('yyyy-MMM')
% ylabel('Sea Surface Temperature [^{\circ}C]','Interpreter','tex')
%% Use pl66 to low-pass filter sst
% with dt=1 (same as default) and T=168 hrs (7 days)
sstf = pl66tn(sst,1,168);
SSTF = timeseries(sstf,datestr(dn)); % Make sstf timeseries
% Plot low-passed filtered data with raw data
figure(2)
plot(SST,'k-')
hold on
plot(SSTF,'--')
title('ERA-5 SST at -35.51, -72.77')
legend('SST_{raw}','Low-passed SST (T=168 hr)')
xlabel('Date')
xtickformat('yyyy-MMM')
ylabel('Sea Surface Temperature [^{\circ}C]','Interpreter','tex')
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

SST0 = timeseries(sst0,datestr(dn)); % Make corresponding timeseries

% Plot climatology
figure(3)
plot(SST0)
title('Climatology from low-passed ERA-5 SST')
xlabel('Date')
xtickformat('yyyy-MMM')
ylabel('Sea Surface Temperature [^{\circ}C]','Interpreter','tex')

%% Take anomaly
sstA = sstf-sst0;
sig = nanstd(sstA);
mean = nanmean(sstA);
SSTA = timeseries(sstA,datestr(dn)); % Make corresponding timeseries

% Plot SST'
figure(4)
yline(mean,'--k')
hold on
yline(mean+2*sig,'--r')
yline(mean-2*sig,'--r')
plot(SSTA)
title("SST' at at -35.51, -72.77")
xlabel('Date')
xtickformat('yyyy-MMM')
ylabel('Sea Surface Temperature [^{\circ}C]','Interpreter','tex')
legend("SST'","\mu_{SST'}","\pm 2 \sigma_{SST'}")

%% First-order approximation of dSST'/dt
dsstdt = sstA(2:end)-sstA(1:end-1); 
DSSTDT = timeseries(dsstdt,datestr(dn(1:end-1))); % Make corresponding timeseries

% Plot of SST' and dt in same figure with 2 y-axes
figure(5)

yyaxis left
plot(SSTA)
title("SST' and dSST'/dt at -35.51, -72.77",'Interpreter','tex')
hold on
yline(0,'--')
xlabel('Date')
xtickformat('yyyy-MMM')
ylabel("SST' [^{\circ}C]",'Interpreter','tex')

yyaxis right
plot(DSSTDT)
hold on
yline(0,'--')
ylabel("dSST'/dt [^{\circ}C/hr]",'Interpreter','tex')
xtickformat('yyyy-MMM')
