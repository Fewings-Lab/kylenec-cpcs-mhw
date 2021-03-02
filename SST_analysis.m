% SST_analysis.m 
% Kylene Cooley 
% 29 Jan 2021 
% A script to carry out analysis of SST anomalies from a netcdf-format
% multi-file dataset

% %% Load dataset from netcdf format
% finfo = ncinfo('Data/Constitucion-SST-2018-2020.nc');
% data1 = squeeze(ncread('Data/Constitucion-SST-1979-1991.nc','sst'));
% data2 = squeeze(ncread('Data/Constitucion-SST-1992-2004.nc','sst'));
% data3 = squeeze(ncread('Data/Constitucion-SST-2005-2017.nc','sst'));
% data4 = ncread('Data/Constitucion-SST-2018-2020.nc','sst');
% data4 = squeeze(data4(:,:,1,:));
% lat = ncread('Data/Constitucion-SST-2018-2020.nc','latitude');
% lon = ncread('Data/Constitucion-SST-2018-2020.nc','longitude');
% time1 = ncread('Data/Constitucion-SST-1979-1991.nc','time');
% time2 = ncread('Data/Constitucion-SST-1992-2004.nc','time');
% time3 = ncread('Data/Constitucion-SST-2005-2017.nc','time');
% time4 = ncread('Data/Constitucion-SST-2018-2020.nc','time');
% 
% % Combine like variables from different files
% sst = cat(1,data1,data2,data3,data4);
% time = cat(1,time1,time2,time3,time4);
% 
% % get datenum vector
% tf = (length(time)/24)+1;
% dn = datenum(1979,1,1:1/24:tf).';
% dn = dn(1:end-1);
% 
% % Save variables as .mat file
% save('ptLavapie_SST.mat','lat','lon','sst','time','dn')

%% Load data from .mat file
clear
load('ptLavapie_SST.mat')
sst = sst-273.15; % Convert sst to degC
ds = datestr(dn);
t = datetime(ds);

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
    mu = mean(sst(values),'omitnan');
    foo(values) = mu;
end 
% Stitch 5 instances of foo together to filter
foo2 = cat(1,foo,foo,foo,foo,foo);

% Filter the working variable
foo3 = pl66tn(foo2,1,240); % apply 10-day filter
% foo4 = pl66tn(foo3,1,hours(years(0.5))); % 6-month (half of a year) low-pass filter
% % take high-pass part of signal to bandpass
% foo5 = foo3-foo4;

% % If we did the high-pass filter first:
% foo3 = pl66tn(foo2,1,hours(years(0.5))); % 6-month (half of a year) low-pass filter
% % take high-pass part of signal to bandpass
% foo4 = foo2-foo3;
% foo5 = pl66tn(foo4,1,240); % apply 10-day filter

% take only the middle portion of filtered climatology
sst0 = foo3(2*length(yd)+1:3*length(yd)); % low-pass filtered only
% sst0 = foo5(2*length(yd)+1:3*length(yd)); % bandpass filtered

% Plot climatology
figure(8)
plot(t,sst0)
% plot(dn(1:(365*24*2)),sst0(1:(365*24*2)))
% Title for low-pass climatology
title('Low-pass filtered climatology of ERA-5 SST')
% % Title for bandpass climatology
% title('Bandpass filtered climatology of ERA-5 SST')
xlabel('Date')
% xlim([dn(1) dn((365*24*2))])
datetick('x','mmmm','keeplimits','keepticks')
ylabel('Sea Surface Temperature [^{\circ}C]','Interpreter','tex')

%% Take anomaly
sstA = sstf-sst0; % using low-pass filtered signal
% sstA = sstb-sst0; % using bandpass filtered signal
% Bandpass SST'
% 10-day low-pass filter
sstA = pl66tn(sstA,1,240); 
% 6-month (half of a year) low-pass filter
hrs = hours(years(0.5)); % define the cutoff frequency in hours
foo6 = pl66tn(sstA,1,hrs); % Evaluate the low-pass filtered signal

% Take high-pass part of signal
sstA = sstA-foo6;
% Replace one window-length with NaNs on each end
sstA(1:2*round(hrs))=NaN;
sstA(end-2*round(hrs):end)=NaN;


sig = nanstd(sstA);
mean = nanmean(sstA);

% Plot SST'
figure(9)
plot(t,sstA)
% % Title for low-pass SST' signal
% title("10-day Low-pass Filtered Sea Surface Temperature Anomaly at 35.51^{\circ}S, 72.77^{\circ}W",'Interpreter','tex')
% Title for bandpass SST' signal
title("10-day to 6-month Bandpass Filtered Sea Surface Temperature Anomaly at 35.51^{\circ}S, 72.77^{\circ}W",'Interpreter','tex')
xlabel('Date')
datetick('x','yyyy-mmm','keeplimits','keepticks')
hold on
yline(mean+2*sig,'--','Color',[0.00,0.45,0.74])
yline(mean-2*sig,'--','Color',[0.00,0.45,0.74])
yline(mean,'-k')
ylabel("SST' [^{\circ}C]",'Interpreter','tex')
legend("SST'","\pm 2 \sigma_{SST'}")

%% First-order approximation of dSST'/dt
dsstdt = sstA(2:end)-sstA(1:end-1); % Using a difference approximation
dsstdtf = pl66tn(dsstdt,1,(3*7*24)); % low-pass filter dSST'/dt with 7-day cutoff
dsstdtf = dsstdtf*24; % Convert to change to degC/day from degC/hr

% Plot of SST' and dt in same figure with 2 y-axes
figure(10)

% Left axis with SST'
yyaxis left
plot(t,sstA)
% % Title for low-pass SST' signal
% title("$\textsf{10-day Low-pass Filtered SST' and 3-week Low-pass Filtered }\frac{\partial \textsf{SST'}}{\partial \textsf{t}}\textsf{ at 35.51}^{\circ}\textsf{S, 72.77}^{\circ}\textsf{W}$",'Interpreter','latex')
% Title for bandpass SST' signal
title("$\textsf{10-day to 6-month Bandpass Filtered SST' and 3-week Low-pass Filtered }\frac{\partial \textsf{SST'}}{\partial \textsf{t}}\textsf{ at 35.51}^{\circ}\textsf{S, 72.77}^{\circ}\textsf{W}$",'Interpreter','latex')
hold on
yline(0,'-k')
xlabel('Date')
datetick('x','yyyy-mmm','keeplimits','keepticks')
yline(mean+2*sig,'--','Color',[0.00,0.45,0.74])
yline(mean-2*sig,'--','Color',[0.00,0.45,0.74])
ylabel("SST' [^{\circ}C]",'Interpreter','tex')

% Right axis with dSST'/dt
yyaxis right
plot(t(1:end-1),dsstdtf)
hold on
yline(0,'-k')
ylabel("$\frac{\partial \textsf{SST'}}{\partial \textsf{t}} \textsf{[}^{\circ}\textsf{C/day]}$",'Interpreter','latex')
datetick('x','yyyy-mmm','keeplimits','keepticks')

% % Limits for IOVWST poster figures
% xlim(datetime(["2015-Jan-01" "2017-Jan-01"]))
% datetick('x','yyyy-mmm','keeplimits')


%% Using matlab functions to find times of all SST' and dSST'dt peaks we are interested in

% Evaluate islocalmax() to find all local maxima of SST' and dSST'dt
% We know there will be at least 1 local max in dSST'dt that leads each
% peak in SST' 

% First, we find all local maxima in each array
SSTlm = islocalmax(sstA);
dSSTlm = islocalmax(dsstdtf);

% Now we will further limit these with SST' that are only above two stds
warmest  = sstA > (mean+2*sig);

% Then select only these values to plot with symbols using logical arrays
% and "AND" operator
sstM = sstA(SSTlm & warmest);
tTm = t(SSTlm & warmest); % times of the local maxima in SST'
dsstM = dsstdtf(dSSTlm);
tdTm  = t(dSSTlm); % times of the local maxima in dSST'

% Selecting events that occur in Dec-Feb
monthnum = dv(SSTlm & warmest, 2); % Pick month number of all of these events
summerSSTa = sstM(monthnum <3 | monthnum==12); % make a vector of these events
summerDates = tTm(monthnum <3 | monthnum==12); % corresponding times

% Reduce dSST' points by limiting to the last nearest date within 24 days of peak SST' dates
maxdate = tTm;
mindate = tTm - days(24);
tdTind = [];

for i = 1:length(maxdate)
    k = find((tdTm > mindate(i) & tdTm < maxdate(i)),1,'last');
    tdTind = cat(1,tdTind,k);
end
% tdTind = ismembertol(datenum(tdTm),datenum(tTm),datenum(days(8))); 
dsstM = dsstM(tdTind);
tdTm = tdTm(tdTind);


% Plot of SST' and dt in same figure with 2 y-axes
figure(11)

yyaxis left
plot(t,sstA,tTm,sstM,'.','MarkerSize',12)
% % Title for low-pass SST' signal
% title("$\textsf{10-day Low-pass Filtered SST' and 3-week Low-pass Filtered }\frac{\partial \textsf{SST'}}{\partial \textsf{t}}\textsf{ at 35.51}^{\circ}\textsf{S, 72.77}^{\circ}\textsf{W}$",'Interpreter','latex')
% Title for bandpass SST' signal
title("$\textsf{10-day to 6-month Bandpass Filtered SST' and 3-week Low-pass Filtered }\frac{\partial \textsf{SST'}}{\partial \textsf{t}}\textsf{ at 35.51}^{\circ}\textsf{S, 72.77}^{\circ}\textsf{W}$",'Interpreter','latex')
hold on
plot(summerDates,summerSSTa,'*','MarkerSize',12,'LineWidth',1)
yline(0,'-k')
xlabel('Date')
datetick('x','yyyy-mmm','keeplimits','keepticks')
yline(mean+2*sig,'--','Color',[0.00,0.45,0.74])
yline(mean-2*sig,'--','Color',[0.00,0.45,0.74])
ylabel("SST' [^{\circ}C]",'Interpreter','tex')

yyaxis right
plot(t(1:end-1),dsstdtf,tdTm,dsstM,'.','MarkerSize',12)
hold on
yline(0,'-k')
ylabel("$\frac{\partial \textsf{SST'}}{\partial \textsf{t}} \textsf{[}^{\circ}\textsf{C/day]}$",'Interpreter','latex')
datetick('x','yyyy-mmm','keeplimits','keepticks')

% % Routines for making decade plots -- repeating these manually
% xlim(datetime(["2010-Jan-01" "2020-Sep-18"]))
% datetick('x','yyyy-mmm','keeplimits')
% title("$\textsf{10-day Low-pass Filtered SST' and 3-week Low-pass Filtered }\frac{\partial \textsf{SST'}}{\partial \textsf{t}}\textsf{ at 35.51}^{\circ}\textsf{S, 72.77}^{\circ}\textsf{W 2010-2020}$",'Interpreter','latex')

