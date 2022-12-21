% upwelling_offshore_ts_comparison.m
% Kylene Cooley
% 28 Apr 2021
% Compares events and peak rates of warming from three separate time
% series: The original "yellow dot" time series near Punta Lavapie
% upwelling center at the highest SST' in the event from the proposal; The
% spatial average from the 1deg sq box 0.5 degrees west of the yellow
% point, which is representative of the upwelling zone; The spatial average
% from a 1deg sq box in the middle of the offshore region of high SST' that
% appears during most "yellow dot" events. 

% 17 June 2021: Change color of orange time series for manuscript Figure 3
% to red. Also changing symbol for times of peak anomalous warming in
% timeseries to filled triangles to differentiate graphically. -- KC

load("sstSwA.mat")

% get indices for lat and lon of either point or center of box
ptLat = [-35.5,-36,-29];        % deg N of dot or center of box
ptLon = [-72.75,-73.75,-77.75]; % deg E of dot or center of box
ind = NaN(2,3);
ind(1,:) = find(ismembertol(lat,ptLat));
ind(2,:) = find(ismembertol(lon,ptLon));
[x1,x2,x3] = matsplit(flip(ind(2,:)));          % sstSwA dim order is lon,lat,time
[y1,y2,y3] = matsplit(circshift(ind(1,:),2));

% time series at yellow dot
ts1 = squeeze(sstSwA(x1,y1,:));

% time series at green box
ts2 = squeeze(mean(sstSwA(x2-2:x2+2,y2-2:y2+2,:),[1 2],'omitnan'));

% time series at offshore box
ts3 = squeeze(mean(sstSwA(x3-2:x3+2,y3-2:y3+2,:),[1 2],'omitnan'));

% remove swath data from workspace to save RAM
clear sstSwA

% take 1st difference approx to time derivative
dts1 = (ts1(2:end)-ts1(1:end-1))*4;         % multiply by 4 for degC/day
dts2 = (ts2(2:end)-ts2(1:end-1))*4; 
dts3 = (ts3(2:end)-ts3(1:end-1))*4; 

% ID events
[sumSST1,sumDates1,sumdSST1,sumdDates1] = SSTevents(time1,ts1,dts1);
[sumSST2,sumDates2,sumdSST2,sumdDates2] = SSTevents(time1,ts2,dts2);
[sumSST3,sumDates3,sumdSST3,sumdDates3] = SSTevents(time1,ts3,dts3);

% plot 3 time series of each variable in same figure with 3 different line
% styles
figure(1) % SST' timeseries
plot(time1,ts1,'-k',sumDates1,sumSST1,'*k','MarkerSize',12,'LineWidth',1)
hold on
plot(time1,ts2,'--b',sumDates2,sumSST2,'*b','MarkerSize',12,'LineWidth',1)
plot(time1,ts3,'-.r',sumDates3,sumSST3,'*r','MarkerSize',12,'LineWidth',1)
yline(0,'-k')
xlabel('Date')
datetick('x','yyyy-mmm','keeplimits','keepticks')
ylabel("SST' [^{\circ}C]",'Interpreter','tex')
legend('Nearshore point','','Spatial mean within nearshore box','','Spatial mean within offshore box')
title("10-day to 6-month bandpass filtered SST'",'Interpreter','latex')

figure(2) % dSST'/dt timeseries
plot(time1(1:end-1),dts1,"-k",sumdDates1,sumdSST1,'*k','MarkerSize',12,'LineWidth',1)
hold on
plot(time1(1:end-1),dts2,"--b",sumdDates2,sumdSST2,'*b','MarkerSize',12,'LineWidth',1)
plot(time1(1:end-1),dts3,"-.r",sumdDates3,sumdSST3,'*r','MarkerSize',12,'LineWidth',1)
yline(0,'-k')
xlabel('Date')
ylabel("$\frac{\partial \textsf{SST'}}{\partial \textsf{t}} \textsf{[}^{\circ}\textsf{C/day]}$",'Interpreter','latex')
datetick('x','yyyy-mmm','keeplimits','keepticks')
legend('Nearshore point','','Spatial mean within nearshore box','','Spatial mean within offshore box')
title("10-day to 6-month bandpass filtered $\frac{\partial \textsf{SST'}}{\partial \textsf{t}}$",'Interpreter','latex')

%% Plotting each location on a separate subplot
figure(2)

ax1 = subplot(3,1,1);
yyaxis left                 % SST' time series
plot(time1,ts1,'-',sumDates1,sumSST1,'*','MarkerSize',12,'LineWidth',1)
hold on
yline(0,'-k')
ylabel("SST' [^{\circ}C]",'Interpreter','tex')
yyaxis right                % dSST'/dt time series
plot(time1(1:end-1),dts1,"-",sumdDates1,sumdSST1,'*','MarkerSize',12,'LineWidth',1)
hold on
yline(0,'-k')
ylabel("$\frac{\partial \textsf{SST'}}{\partial \textsf{t}} \textsf{[}^{\circ}\textsf{C/day]}$",'Interpreter','latex')
datetick('x','yyyy-mmm','keeplimits','keepticks')
title("Nearshore point (N_{events}=31)")

ax2 = subplot(3,1,2);
yyaxis left                 % SST' time series
plot(time1,ts2,'--',sumDates2,sumSST2,'*','MarkerSize',12,'LineWidth',1)
hold on
yline(0,'-k')
ylabel("SST' [^{\circ}C]",'Interpreter','tex')
yyaxis right                % dSST'/dt time series
plot(time1(1:end-1),dts2,"--",sumdDates2,sumdSST2,'*','MarkerSize',12,'LineWidth',1)
hold on
yline(0,'-k')
ylabel("$\frac{\partial \textsf{SST'}}{\partial \textsf{t}} \textsf{[}^{\circ}\textsf{C/day]}$",'Interpreter','latex')
datetick('x','yyyy-mmm','keeplimits','keepticks')
title("Spatial mean within nearshore box (N_{events}=36)")

ax3 = subplot(3,1,3);
yyaxis left                 % SST' time series
plot(time1,ts3,'-.',sumDates3,sumSST3,'*','MarkerSize',12,'LineWidth',1)
hold on
yline(0,'-k')
xlabel('Date')
datetick('x','yyyy-mmm','keeplimits','keepticks')
ylabel("SST' [^{\circ}C]",'Interpreter','tex')
yyaxis right                % dSST'/dt time series
plot(time1(1:end-1),dts3,"-.",sumdDates3,sumdSST3,'*','MarkerSize',12,'LineWidth',1)
hold on
yline(0,'-k')
ylabel("$\frac{\partial \textsf{SST'}}{\partial \textsf{t}} \textsf{[}^{\circ}\textsf{C/day]}$",'Interpreter','latex')
datetick('x','yyyy-mmm','keeplimits','keepticks')
title("Spatial mean within offshore box (N_{events}=41)")

sgtitle("10-day to 6-month bandpass filtered SST' and $\frac{\partial \textsf{SST'}}{\partial \textsf{t}}$",'Interpreter','latex')
linkaxes([ax1 ax2 ax3],'x')

%% simple timeseries for nearshore box mean
sig_box = std(ts2,'omitnan');

figure()
plot(time1,ts2,'b',sumDates2,sumSST2,'*b','MarkerSize',12,'LineWidth',1)
hold on 
yline(0,'-k')
yline(2*sig_box,'--b')
yline(-2*sig_box,'--b')
xlabel('Date')
datetick('x','yyyy-mmm','keeplimits','keepticks')
ylabel("SST' [^{\circ}C]",'Interpreter','tex')
legend(["SST'";'';'';'\pm 2\sigma limits'])
title("10-day to 6-month bandpass filtered SST'",'Interpreter','latex')

%% time series without stars
sig_box = std(ts2,'omitnan');

figure()
plot(time1,ts2,'b')
hold on 
yline(0,'-k')
yline(2*sig_box,'--b')
yline(-2*sig_box,'--b')
xlabel('Date')
datetick('x','yyyy-mmm','keeplimits','keepticks')
ylabel("SST' [^{\circ}C]",'Interpreter','tex')
legend(["SST'";'';'\pm 2\sigma limits'])
title("10-day to 6-month bandpass filtered SST'",'Interpreter','latex')
%% simple dSST'/dt for nearshore box
h1 = figure();
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')

nexttile
%ax1 = subplot(2,1,1);
colororder({'#0072BD','r'})
yyaxis left      % SST' time series
plot(time1,ts2,sumDates2,sumSST2,'*','MarkerSize',8,'LineWidth',0.7)
hold on
yline(0,'-k')
ylim([-2.5 2.5])
ylabel("SST' [^{\circ}C]",'Interpreter','tex','FontSize',8)
yyaxis right                % dSST'/dt time series
plot(time1(1:end-1),dts2,"--r",sumdDates2,sumdSST2,'^','MarkerFaceColor','r','MarkerSize',6,'LineWidth',0.7)
hold on
yline(0,'-k')
ylim([-0.5 0.5])
ylabel("$\partial \textsf{SST'}/\partial \textsf{t} \textsf{ [}^{\circ}\textsf{C day}^{-1}\textsf{]}$",'Interpreter','latex','FontSize',8)
datetick('x','yyyy-mmm','keeplimits','keepticks')
%title("10-day to 6-month bandpass filtered SST' and $\frac{\partial SST'}{\partial t}$",'Interpreter','latex')
xlabel('Date','FontSize',8)
set(gca,'FontSize',8)

nexttile
% ax2 = subplot(2,1,2);
colororder({'#0072BD','r'})
yyaxis left      % SST' time series
plot(time1,ts2,sumDates2,sumSST2,'*','MarkerSize',8,'LineWidth',0.7)
hold on
yline(0,'-k')
ylim([-2.5 2.5])
ylabel("SST' [^{\circ}C]",'Interpreter','tex','FontSize',8)
yyaxis right                % dSST'/dt time series
plot(time1(1:end-1),dts2,"--",sumdDates2,sumdSST2,'^','MarkerSize',6,'MarkerFaceColor','r','LineWidth',0.7)
hold on
yline(0,'-k')
ylim([-0.5 0.5])
ylabel("$\partial \textsf{SST'}/\partial \textsf{t} \textsf{ [}^{\circ}\textsf{C day}^{-1}\textsf{]}$",'Interpreter','latex','FontSize',8)
xlim([datetime('2007-11-01') datetime('2010-03-01')])
datetick('x','yyyy-mmm','keeplimits','keepticks')
%title("10-day to 6-month bandpass filtered SST' and $\frac{\partial SST'}{\partial t}$",'Interpreter','latex')
xlabel('Date','FontSize',8)
set(gca,'FontSize',8)

h1.Units = 'centimeters';
h1.Position = [0 0 19 12];
h1.PaperUnits = 'centimeters';
h1.PaperPosition = [0 0 19 12];

%exportgraphics(h1,"SST'-dSST'-box-avg-simple-v6.pdf")