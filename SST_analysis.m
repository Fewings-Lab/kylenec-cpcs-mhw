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
t = datetime(ds);
% SST = timeseries(sst,ds); % Make sst timeseries object
% % Plot raw data
% figure(1)
% plot(t,sst)
% title('ERA-5 Sea Surface Temperature at 35.51^{\circ}S, 72.77^{\circ}W')
% xlabel('Date')
% datetick('x','yyyy-mmm','keeplimits','keepticks')
% ylabel('SST [^{\circ}C]','Interpreter','tex')
%% Use pl66 to low-pass filter sst
% with dt=1 (same as default) and T=168 hrs (7 days)
sstf = pl66tn(sst,1,240);
% SSTF = timeseries(sstf,ds); % Make sstf timeseries
% Plot low-passed filtered data with raw data
figure(7)
plot(t,sst,'Color',[0.8,0.8,0.8],'LineWidth',2)
hold on
plot(t,sstf,'k','LineWidth',1)
title('ERA-5 Sea Surface Temperature at 35.51^{\circ}S, 72.77^{\circ}W','Interpreter','tex')
legend('SST_{raw}','10-day low-pass filtered SST')
xlabel('Date')
datetick('x','yyyy-mmm','keeplimits','keepticks')
ylabel('SST [^{\circ}C]','Interpreter','tex')

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
foo3 = pl66tn(foo2,1,240); % apply 7-day filter

% foo4 = pl66tn(foo3,1,(30*24)); % 30-day filter
% % take high-pass part of signal to bandpass
% foo5 = foo3-foo4;

sst0 = foo3(length(yd)+1:2*length(yd)); % take only the middle portion of filtered climatology

% SST0 = timeseries(sst0,ds); % Make corresponding timeseries

% Plot climatology
figure(8)
plot(t,sst0)
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
% SSTA = timeseries(sstA,ds); % Make corresponding timeseries

% Plot SST'
figure(9)
plot(t,sstA)
title("10-day Low-pass Filtered Sea Surface Temperature Anomaly at 35.51^{\circ}S, 72.77^{\circ}W",'Interpreter','tex')
xlabel('Date')
datetick('x','yyyy-mmm','keeplimits','keepticks')
hold on
yline(mean+2*sig,'--','Color',[0.00,0.45,0.74])
yline(mean-2*sig,'--','Color',[0.00,0.45,0.74])
yline(mean,'-k')
ylabel("SST' [^{\circ}C]",'Interpreter','tex')
legend("SST'","\pm 2 \sigma_{SST'}")

%% First-order approximation of dSST'/dt
dsstdt = sstA(2:end)-sstA(1:end-1); 
dsstdtf = pl66tn(dsstdt,1,(3*7*24));
dsstdtf = dsstdtf*24; % Convert to change to degC/day from degC/hr
% DSSTDT = timeseries(dsstdtf,ds(1:end-1,:)); % Make corresponding timeseries

% Plot of SST' and dt in same figure with 2 y-axes
figure(10)

yyaxis left
plot(t,sstA)
title("$\textsf{10-day Low-pass Filtered SST' and }\frac{\partial \textsf{SST'}}{\partial \textsf{t}}\textsf{ at 35.51}^{\circ}\textsf{S, 72.77}^{\circ}\textsf{W}$",'Interpreter','latex')
hold on
yline(0,'-k')
xlabel('Date')
datetick('x','yyyy-mmm','keeplimits','keepticks')
yline(mean+2*sig,'--','Color',[0.00,0.45,0.74])
yline(mean-2*sig,'--','Color',[0.00,0.45,0.74])
ylabel("SST' [^{\circ}C]",'Interpreter','tex')

yyaxis right
plot(t(1:end-1),dsstdtf)
hold on
yline(0,'-k')
ylabel("$\frac{\partial \textsf{SST'}}{\partial \textsf{t}} \textsf{[}^{\circ}\textsf{C/day]}$",'Interpreter','latex')
datetick('x','yyyy-mmm','keeplimits','keepticks')

% Limits for IOVWST poster figures
xlim(datetime(["2015-Jan-01" "2017-Jan-01"]))
datetick('x','yyyy-mmm','keeplimits')