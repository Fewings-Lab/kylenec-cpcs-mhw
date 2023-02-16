% explore_winds.m
% A script to store all the snippets of code that made plots that I have
% but ultimately turned out to have the wrong preprocessing applied.
% Kylene Cooley
% 25 Aug 2021


%% Anomaly of wind stress magnitude and northward component only
windmagA = windstrmag - windclim;
vwindA = y - vclim;
uwindA = x - uclim;

figure(1)
plot(time(14000:14365),squeeze(windclim(100,30,14000:14365)),'k-',time(14000:14365),squeeze(vclim(100,30,14000:14365)),'r-.',time(14000:14365),squeeze(uclim(100,30,14000:14365)),'b:')
legend('wind stress magnitude','meridional wind stress','zonal wind stress')
datetick('keeplimits')
xlabel('Date')
ylabel('Turbulent wind stress [N/m^2]')
title('Surface wind stress climatology at point 25.5^\circ S, 82.75^\circ W')
% % apply a bandpass filter -- since we applied a bandpass filter to the time
% series, we may not need this.
% [windmagA,~,~] = bandpass(windmagA,24,240,0.5);
%%
figure(2)
plot(time(14000:14365),squeeze(windmagA(100,30,14000:14365)),'k-',time(14000:14365),squeeze(vwindA(100,30,14000:14365)),'r-.',time(14000:14365),squeeze(uwindA(100,30,14000:14365)),'b:')
legend('wind stress magnitude','meridional wind stress','zonal wind stress')
datetick('keeplimits')
xlabel('Date')
ylabel('Turbulent wind stress [N/m^2]')
title('Surface wind stress anomaly at point 25.5^\circ S, 82.75^\circ W')

%%
figure(3)
plot(time(14000:14365),squeeze(windstrmag(100,30,14000:14365)),'k-',time(14000:14365),squeeze(y(100,30,14000:14365)),'r-.',time(14000:14365),squeeze(x(100,30,14000:14365)),'b:')
legend('wind stress magnitude','meridional wind stress','zonal wind stress')
datetick('keeplimits')
xlabel('Date')
ylabel('Turbulent wind stress [N/m^2]')
title('Surface wind stress at point 25.5^\circ S, 82.75^\circ W')

%%
latidx = find(lat==-35);
lonidx = find(lon==(360-75));

figure(4)
plot(time(14000:14365),squeeze(windclim(latidx,lonidx,14000:14365)),'k-',time(14000:14365),squeeze(vclim(latidx,lonidx,14000:14365)),'r-.',time(14000:14365),squeeze(uclim(latidx,lonidx,14000:14365)),'b:')
legend('wind stress magnitude','meridional wind stress','zonal wind stress')
datetick('keeplimits')
xlabel('Date')
ylabel('Turbulent wind stress [N/m^2]')
title('Surface wind stress climatology at point 35^\circ S, 75^\circ W')

figure(5)
plot(time(14000:14365),squeeze(windmagA(latidx,lonidx,14000:14365)),'k-',time(14000:14365),squeeze(vwindA(latidx,lonidx,14000:14365)),'r-.',time(14000:14365),squeeze(uwindA(latidx,lonidx,14000:14365)),'b:')
legend('wind stress magnitude','meridional wind stress','zonal wind stress')
datetick('keeplimits')
xlabel('Date')
ylabel('Turbulent wind stress [N/m^2]')
title('Surface wind stress anomaly at point 35^\circ S, 75^\circ W')

figure(6)
plot(time(14000:14365),squeeze(windstrmag(latidx,lonidx,14000:14365)),'k-',time(14000:14365),squeeze(y(latidx,lonidx,14000:14365)),'r-.',time(14000:14365),squeeze(x(latidx,lonidx,14000:14365)),'b:')
legend('wind stress magnitude','meridional wind stress','zonal wind stress')
datetick('keeplimits')
xlabel('Date')
ylabel('Turbulent wind stress [N/m^2]')
title('Surface wind stress at point 35^\circ S, 75^\circ W')
%% Save workspace as -.mat file
save('windstress_era5.mat','-append')