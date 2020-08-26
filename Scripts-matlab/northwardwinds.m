% northwardwinds.m
% Author: Kylene Cooley
% Created 3 May 2020
% Edited 29 June 2020: Originally from termproject.m which I wrote to look at high-frequency v. low-frequency(~diurnal) variability in the northward component of 10-m wind velocity. Used for my term project in OC 682 (Data Analysis in Time and Space Domain).

% Read in 10-m wind velocity northward component and lat/lon data
% finfo = ncinfo('V10_2010_2020.nc'); % Used this to find what names of
% variables were
v10 = ncread('V10_2010_2020.nc','v10'); % 10-m wind velocity northward component, 9x9x2x90934 double
lat = ncread('V10_2010_2020.nc','latitude'); % Latitude values for each grid pt., 9x1 double
lon = ncread('V10_2010_2020.nc','longitude'); % Longitude values, 9x1 double
t = double(ncread('V10_2010_2020.nc','time')); % Time vector, 90934x1 double

% Pull out data for the grid point near Pt. Lengua de Vaca:
dat = rmmissing(squeeze(v10(5,7,:))); % Raw time series, 90934x1 double

% Make new date number vector for sorting data later:
% disp(datestr(t(1))) % Used this to check starting and ending date values
% disp(datestr(t(end)))
Datestart = datenum(2010,01,01,00,00,00);
Dateend = datenum(2020,05,16,22,30,00);
dates = datenum(2010, 01, 01,[0:((Dateend-Datestart)*24)-1].',0,0); % Hourly date numbers, 90934x1 double

% Plot of unfiltered time series
figure(1)
plot(dates,dat)
datetick()
hold on
yline(0,'--');
xlabel('Date')
ylabel('V_{10}(t) [m/s]')
title('Northward Component of 10-meter Wind Velocity Near Lengua de Vaca Station')
xlim([dates(1) dates(end)])

% Low pass filter time series to obtain separate high and low frequency signals
datLo = pl66tn(dat);
datHi = dat-datLo;

% Plot these high and low frequency time series
figure(3)
subplot(2,1,1)
plot(dates,datLo)
datetick()
title('Low-Pass Filtered Northward Component of 10-meter Wind Velocity')
ylabel('V_{10}(t) [m/s]')
xlabel('Date')
xlim([dates(1) dates(end)])
subplot(2,1,2)
plot(dates,datHi)
datetick()
title('High-Pass Filtered Northward Component of 10-meter Wind Velocity')
ylabel('V_{10}(t) [m/s]')
xlabel('Date')
xlim([dates(1) dates(end)])

% To sort data by month of measurement, create date number vector
Dates = datevec(dates); % Date number vector, 90934x6 double

% Sorting raw and both filtered sets by month
datmat = [dat, datLo, datHi]; % Concatenating the raw and filtered data for easier loops
monthDatmat = NaN.*ones(3,12,90934); % NaN matrix to hold data separated by month
for i = 1:3
    wheat = datmat(:,i); % Working variable for column of data used in this loop
    for j = 1:12
        % Assigning spaces in the NaN matrix to hold data separated by
        % month and time series
        monthDatmat(i,j,(1:length(wheat(Dates(:,2) == j)))) = wheat(Dates(:,2) == j);
    end
end

% Get mean for each month across rows
mDat = nanmean(monthDatmat(1,:,:),3); % A way to check I did this right
mLo = nanmean(monthDatmat(2,:,:),3);
mHi = nanmean(monthDatmat(3,:,:),3);
% Get standard deviation for each month across rows
stdLo = nanstd(monthDatmat(2,:,:),0,3);
stdHi = nanstd(monthDatmat(3,:,:),0,3);

% Making vectors of 2 cycles to plot
monvec = [1:12,1:12];
monstr = num2str(monvec','%i');
meanvecL = cat(2,mLo,mLo);
stdvecL = cat(2,stdLo,stdLo);
meanvecH = cat(2,mHi,mHi);
stdvecH = cat(2,stdHi,stdHi);

% Get confidence intervals on the mean to add to graph
N = length(dat);
CILo = (1.96/sqrt(10)).*stdvecL;
CIHi = (1.96/sqrt(10)).*stdvecH;

% Graph of monthly means with standard deviation errorbars and confidence intervals
% Lower frequencies
figure(4)
errorbar(meanvecL,stdvecL)
hold on
plot(meanvecL+CILo,'r--','LineWidth',1)
plot(meanvecL-CILo,'r--','LineWidth',1)
title('Low-Pass Filtered Monthly Mean Northward Component of Wind Velocity')
ylabel('V_{10}(t) [m/s]')
xticks(1:24)
xticklabels(monstr)
xlabel('Month Number')
legend('Monthly Mean and Standard Deviation','95% Confidence Interval')

% Higher frequencies
figure(5)
errorbar(meanvecH,stdvecH)
hold on
plot(meanvecH+CIHi,'r--','LineWidth',1)
plot(meanvecH-CIHi,'r--','LineWidth',1)
title('Diurnal and Higher Frequency Monthly Mean Northward Component of Wind Velocity')
ylabel('V_{10}(t) [m/s]')
xticks(1:24)
xticklabels(monstr)
xlabel('Month Number')
legend('Monthly Mean and Standard Deviation','95% Confidence Interval')

% Variability explained by high and low-pass filtered data individually
% Variance for total time series and filtered sets
vary = 1/N*sum((dat-nanmean(dat)).^2);
varyL = 1/N*sum((datLo-nanmean(datLo)).^2);
varyH = 1/N*sum((datHi-nanmean(datHi)).^2);

% Variability explained by high and low pass signals
varbL = varyL/vary;
varbH = varyH/vary;

% Checking amount of cross terms that explain variability
% Variance of cross terms minus low and high filtered variance
% varyDif = 1/length(dat)*sum(2.*(datLo.*datHi-datLo*nanmean(dat)+datLo*nanmean(datLo)-datHi*nanmean(dat)+datHi*nanmean(datHi))+(nanmean(dat))^2-(nanmean(datLo))^2-(nanmean(datHi))^2);
% varbDif = varyDif/vary;

% High and low frequency variability by month
monthVarb = NaN.*ones(3,12); % first column will be changed to the sum of variability from yH and yL after loop; 
monthVary = NaN.*ones(3,12); % Variance matrix for [total monthly signal, low freq, high freq]
for i = 1:3
    for j = 1:12
        % Assigning spaces in the NaN matrix to hold data separated by
        % month and time series
        monthVary(i,j) = 1/N*nansum((monthDatmat(i,j,:)-nanmean(monthDatmat(i,j,:))).^2);
        if i > 1
            monthVarb(i,j) = monthVary(i,j)/monthVary(1,j);
        end
    end
end
monthVarb(1,:) = monthVarb(2,:) + monthVarb(3,:);

% Plot low, high, and combined variability
figure(6)
plot(monthVarb(1,:))
hold on
plot(monthVarb(2,:))
plot(monthVarb(3,:))
legend('Combined high and low frequencies','Low-passed frequencies','High-passed frequencies')
title('Variability by Month')
xlabel('Month number')
ylabel('Fraction variability explained')
