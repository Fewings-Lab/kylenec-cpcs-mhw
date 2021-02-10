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
ds = datestr(dn);
SST = timeseries(sst,ds); % Make sst timeseries object
% % Plot raw data
% figure(1)
% plot(sst)
% title('ERA-5 SST at -35.51, -72.77')
% xlabel('Date')
% datetick('x','yyyy-mmm','keeplimits','keepticks')
% ylabel('Sea Surface Temperature [^{\circ}C]','Interpreter','tex')
%% Use pl66 to low-pass filter sst
% with dt=1 (same as default) and T=168 hrs (7 days)
sstf = pl66tn(sst,1,168);
SSTF = timeseries(sstf,ds); % Make sstf timeseries
% Plot low-passed filtered data with raw data
figure(2)
plot(SST,'k-')
hold on
plot(SSTF,'--')
title('ERA-5 SST at -35.51, -72.77')
legend('SST_{raw}','7-day low-pass filtered SST')
xlabel('Date')
datetick('x','yyyy-mmm','keeplimits','keepticks')
ylabel('Sea Surface Temperature [^{\circ}C]','Interpreter','tex')

% %% Bandpass sst
% % 30 day low pass filter
% sstf2 = pl66tn(sstf,1,(30*24));
% % take high-pass part of signal
% sstb = sstf-sstf2;
% % make timeseries
% SSTB = timeseries(sstb,ds);
%% Obtain climatology
% Convert datenums to date vector [yy mm dd hh mm ss]
dv = datevec(dn); 

% Calculate yearday by subtracting the datenum of Jan 1 of that year 
% (and add 1 so the first day of the year is yearday 1 not 0)
yd = dn - datenum(dv(:,1),1,1) + 1; 

% Vectors to use for matching times and storing climatology
yhr = 1:1/24:(367-1/24);
foo = NaN(length(yd),1);

% loop by values
for i=1:366*24
    values = ismembertol(yd,yhr(i));
    mu = nanmean(sst(values));
    foo(values) = mu;
end 
% Stitch 3 instances of foo together to filter
foo2 = cat(1,foo,foo,foo);
foo3 = pl66tn(foo2,1,168); % apply 7-day filter

% foo4 = pl66tn(foo3,1,(30*24)); % 30-day filter
% % take high-pass part of signal to bandpass
% foo5 = foo3-foo4;

sst0 = foo3(length(yd)+1:2*length(yd)); % take only the middle portion of filtered climatology

SST0 = timeseries(sst0,ds); % Make corresponding timeseries

% Plot climatology
figure(3)
plot(SST0)
% plot(dn(1:(365*24*2)),sst0(1:(365*24*2)))
title('Low-pass filtered climatology of ERA-5 SST')
xlabel('Date')
% xlim([dn(1) dn((365*24*2))])
datetick('x','mmmm','keeplimits','keepticks')
ylabel('Sea Surface Temperature [^{\circ}C]','Interpreter','tex')

%% Take anomaly
sstA = sstf-sst0;
sig = nanstd(sstA);
mean = nanmean(sstA);
SSTA = timeseries(sstA,ds); % Make corresponding timeseries

% Plot SST'
figure(4)
plot(SSTA)
title("SST' at -35.51, -72.77")
xlabel('Date')
datetick('x','yyyy-mmm','keeplimits','keepticks')
hold on
yline(mean,'--k')
yline(mean+2*sig,'--r')
yline(mean-2*sig,'--r')
ylabel('Sea Surface Temperature [^{\circ}C]','Interpreter','tex')
legend("SST'","\mu_{SST'}","\pm 2 \sigma_{SST'}")

%% First-order approximation of dSST'/dt
dsstdt = sstA(2:end)-sstA(1:end-1); 
dsstdtf = pl66tn(dsstdt,1,(3*7*24));
DSSTDT = timeseries(dsstdtf,ds(1:end-1,:)); % Make corresponding timeseries

% Plot of SST' and dt in same figure with 2 y-axes
figure(5)

yyaxis left
plot(SSTA)
title("$\textsf{SST' and }\frac{\partial \textsf{SST'}}{\partial \textsf{t}}\textsf{ at -35.51, -72.77}$",'Interpreter','latex')
hold on
yline(0,'--')
xlabel('Date')
datetick('x','yyyy-mmm','keeplimits','keepticks')
yline(mean+2*sig,'--r')
yline(mean-2*sig,'--r')
ylabel("SST' [^{\circ}C]",'Interpreter','tex')

yyaxis right
plot(DSSTDT)
hold on
yline(0,'--')
ylabel("$\frac{\partial \textsf{SST'}}{\partial \textsf{t}} \textsf{[}^{\circ}\textsf{C/hr]}$",'Interpreter','latex')
datetick('x','yyyy-mmm','keeplimits','keepticks')
