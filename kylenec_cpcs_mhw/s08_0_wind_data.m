% s08_0_wind_data.m
% Kylene Cooley
% 5 Aug 2021

% 1) Initial variables defining our research area and time period
latlim = [-50 -15];
lonlim = 360 + [-90 -70];
stdate = '01/01/1979';
endate = '12/31/2020';

go = isfile('windstress_era5_v3.mat');

if go == 0
    % Collect unfiltered eastward (u) and northward (v) wind stress fields
%     uwindstr = read_era5_new('sfc','daily',stdate,endate,{'ustr_sfc'},lonlim,latlim);   
%     vwindstr = read_era5_new('sfc','daily',stdate,endate,{'vstr_sfc'},lonlim,latlim);
    wstrm_anomaly = read_era5_new('sfc','daily',stdate,endate,{'strm_sfc'},lonlim,latlim,'anomaly',1);
    save('windstress_era5_v3.mat','-v7.3') 
end
     
%% 2) Magnitude of the wind stress in N/m^2
% Load structs created from subset of ERA5 wind stress
load('windstress_era5_v3.mat')

% give useful variables short names
% x = uwindstr.ustr_sfc;
% y = vwindstr.vstr_sfc;
time = wstrm_anomaly.time;
lat = wstrm_anomaly.lat;
lon = wstrm_anomaly.lon;

% use the norm to get magnitude
% windstrmag = sqrt((x.^2)+(y.^2));
windstrmag = wstrm_anomaly.anomaly_strm_sfc;
% %% 3) Wind stress magnitude in SST' detection (nearshore) box
% [foo2_mag, foo2_u, foo2_v] = nearshoreboxmean(lat,lon,windstrmag,x,y);
% 
% figure()
% ax1 = subplot(3,1,1);
% plot(time,foo2_mag,'k')
% title('Wind stress magnitude')
% datetick('keeplimits')
% ax2 = subplot(3,1,2);
% plot(time,foo2_u,'b--')
% title('Zonal wind stress component')
% datetick('keeplimits')
% ylabel('Wind stress [N/m^2]')
% ax3 = subplot(3,1,3);
% plot(time,foo2_v,'r-.')
% datetick('keeplimits')
% xlabel('Date')
% title('Meridional wind stress component')
% sgtitle('Mean wind stress in area bounded by 35.5^\circS, 73.25^\circW and 36.5^\circS, 74.25^\circW')
% linkaxes([ax1 ax2 ax3],'x')
% 
% % saveas(gcf,'ts_mean_wind_stress_raw.png')
% 
% %% 4) Climatology of wind stress mag, u , v without a filter and plot
% [windclim,~] = clim1y3d(windstrmag,time,24);
% [vclim,~] = clim1y3d(y,time,24);
% [uclim,~] = clim1y3d(x,time,24);
% 
% [clim_mag,clim_u,clim_v] = nearshoreboxmean(lat,lon,windclim,uclim,vclim);
% 
% figure()
% ax1 = subplot(3,1,1);
% plot(time(13600:13965),clim_mag(13600:13965),'k')
% title('Wind stress magnitude')
% datetick('keeplimits')
% ax2 = subplot(3,1,2);
% plot(time(13600:13965),clim_u(13600:13965),'b--')
% title('Zonal wind stress component')
% datetick('keeplimits')
% ylabel('Wind stress [N/m^2]')
% ax3 = subplot(3,1,3);
% plot(time(13600:13965),clim_v(13600:13965),'r-.')
% datetick('keeplimits')
% xlabel('Date')
% title('Meridional wind stress component')
% sgtitle('Mean wind stress climatology in area bounded by 35.5^\circS, 73.25^\circW and 36.5^\circS, 74.25^\circW')
% linkaxes([ax1 ax2 ax3],'x')
% 
% % saveas(gcf,'ts_clim_wind_stress_raw.png')
% %% 5) Low-pass (Tc = 10 days) filtered climatology and plot
% cmLP = pl66tn3D(windclim,24,240);
% cuLP = pl66tn3D(uclim,24,240);
% cvLP = pl66tn3D(vclim,24,240);
% 
% [cmLP_box,cuLP_box,cvLP_box] = nearshoreboxmean(lat,lon,cmLP,cuLP,cvLP);
% 
% figure()
% ax1 = subplot(3,1,1);
% plot(time(13600:13965),cmLP_box(13600:13965),'k')
% title('Wind stress magnitude')
% datetick('keeplimits')
% ax2 = subplot(3,1,2);
% plot(time(13600:13965),cuLP_box(13600:13965),'b--')
% title('Zonal wind stress component')
% datetick('keeplimits')
% ylabel('Wind stress [N/m^2]')
% ax3 = subplot(3,1,3);
% plot(time(13600:13965),cvLP_box(13600:13965),'r-.')
% datetick('keeplimits')
% xlabel('Date')
% title('Meridional wind stress component')
% sgtitle('Mean wind stress climatology low-pass filtered in area bounded by 35.5^\circS, 73.25^\circW and 36.5^\circS, 74.25^\circW')
% linkaxes([ax1 ax2 ax3],'x')
% 
% % saveas(gcf,'ts_clim_wind_stress_10dLP.png')
% %% 6) Band-pass filtered (Tc,low = 10d, Tc,high = 6mo) climatology and plot
% % Take 6 mo. low pass filter of climatology
% cmLPLP = pl66tn3D(cmLP,24,hours(years(0.5)));
% cuLPLP = pl66tn3D(cuLP,24,hours(years(0.5)));
% cvLPLP = pl66tn3D(cvLP,24,hours(years(0.5)));
% 
% % Take difference of LPLP from LP to get the bandpass filter
% cmBP = cmLP - cmLPLP;
% cuBP = cuLP - cuLPLP;
% cvBP = cvLP - cvLPLP;
% 
% % We will plot just the average in the nearshore box
% [cmBP_box,cuBP_box,cvBP_box] = nearshoreboxmean(lat,lon,cmBP,cuBP,cvBP);
% 
% figure()
% ax1 = subplot(3,1,1);
% plot(time(13600:13965),cmBP_box(13600:13965),'k')
% title('Wind stress magnitude')
% datetick('keeplimits')
% ax2 = subplot(3,1,2);
% plot(time(13600:13965),cuBP_box(13600:13965),'b--')
% title('Zonal wind stress component')
% datetick('keeplimits')
% ylabel('Wind stress [N/m^2]')
% ax3 = subplot(3,1,3);
% plot(time(13600:13965),cvBP_box(13600:13965),'r-.')
% datetick('keeplimits')
% xlabel('Date')
% title('Meridional wind stress component')
% sgtitle('Mean wind stress climatology bandpass filtered in area bounded by 35.5^\circS, 73.25^\circW and 36.5^\circS, 74.25^\circW')
% linkaxes([ax1 ax2 ax3],'x')
% 
% % saveas(gcf,'ts_clim_wind_stress_10d6moBP.png')
%% 7) Plotting climatology of both filtered and unfiltered on top of wind stress magnitude
% Make a bandpass filtered wind stress magnitude
magLP = pl66tn3D(windstrmag,24,240);
magLPLP = pl66tn3D(magLP,24,hours(years(0.5)));
magBP = magLP - magLPLP;

% Mean bandpass filtered wind stress magnitude in the nearshore box
% [magBP_box,mag_box,~] = nearshoreboxmean(lat,lon,magBP,windstrmag,cvBP);

% %% 8) I will show both in the same figure with 2x1 subplots
% figure()
% ax1 = subplot(2,1,1);
% plot(time(13600:13965),mag_box(13600:13965),'k',time(13600:13965),clim_mag(13600:13965),'r-.')
% legend('Wind stress magnitude','Climatology')
% title('Unfiltered wind stress')
% datetick('keeplimits')
% ylabel('Wind stress [N/m^2]')
% ax2 = subplot(2,1,2);
% plot(time(13600:13965),magBP_box(13600:13965),'k',time(13600:13965),cmBP_box(13600:13965),'r-.')
% legend('Wind stress magnitude','Climatology')
% title('Bandpass filtered wind stress')
% datetick('keeplimits')
% ylabel('Wind stress [N/m^2]')
% xlabel('Date')
% sgtitle('Mean wind stress in area bounded by 35.5^\circS, 73.25^\circW and 36.5^\circS, 74.25^\circW')
% linkaxes([ax1 ax2],'x')
% 
% % saveas(gcf,'ts_clim_wind_stress_magBP.png')
% %% 9) plot filtered and unfiltered magnitude on same axis (then same with climatology)
% 
% figure()
% ax1 = subplot(2,1,1);
% plot(time(13600:13965),mag_box(13600:13965),'k',time(13600:13965),magBP_box(13600:13965),'r-.')
% legend('Unfiltered wind stress','10-day to 6-month bandpass filtered ')
% title('Wind stress magnitude')
% datetick('keeplimits')
% ylabel('Wind stress [N/m^2]')
% ax2 = subplot(2,1,2);
% plot(time(13600:13965),clim_mag(13600:13965),'k',time(13600:13965),cmBP_box(13600:13965),'r-.')
% legend('Unfiltered wind stress','10-day to 6-month bandpass filtered ')
% title('Climatology of wind stress magnitude')
% datetick('keeplimits')
% ylabel('Wind stress [N/m^2]')
% xlabel('Date')
% sgtitle('Mean wind stress in area bounded by 35.5^\circS, 73.25^\circW and 36.5^\circS, 74.25^\circW')
% linkaxes([ax1 ax2],'x')
% 
% % saveas(gcf,'ts_wind_stress_climBP_and_magBP.png')
% %% 10) BP filter applied to climatology vs. climatology from BP-filtered magnitude
% [cBPm,~] = clim1y3d(magBP,time,24);
% 
% % Mean bandpass filtered wind stress magnitude in the nearshore box
% [cBPm_box,~,~] = nearshoreboxmean(lat,lon,cBPm,windstrmag,cvBP);
% %% 11) plot BP clim v. clim of BP mag
% figure()
% plot(time(13600:13965),cmBP_box(13600:13965),'b--',time(13600:13965),cBPm_box(13600:13965),'r-.')
% legend('Bandpass-filtered climatology','Climatology from bandpass-filtered magnitude')
% title('Climatology of wind stress magnitude')
% datetick('keeplimits')
% ylabel('Wind stress [N/m^2]')
% xlabel('Date')
% 
% % saveas(gcf,'ts_wind_stress_climBP_compare.png')
%% 12) Save .mat files for each of the wind stress quantities
% save('windstrmag.mat','cBPm','cBPm_box','clim_mag','cmBP','cmBP_box','cmLP','cmLP_box','cmLPLP','mag_box','magBP','magBP_box','magLP','magLPLP','windstrmag','time','lat','lon','windclim','-v7.3')
% save('windstru.mat','clim_u','cuBP','cuBP_box','cuLP','cuLP_box','cuLPLP','uclim','uwindA','x','time','lat','lon','-v7.3')
% save('windstrv.mat','clim_v','cvBP','cvBP_box','cvLP','cvLP_box','cvLPLP','vclim','vwindA','y','time','lat','lon','-v7.3')
save('windstrmag2.mat','magBP','magLP','windstrmag','time','lat','lon','-v7.3')

%% 13) 10d to 6m bandpass filtered wind stress magnitude anomaly: composite event map

load windstrmag2.mat
load 0_05_level_lines.mat
load 0_035_residual_lines.mat

cmBP = zeros(size(magBP));
magBP_cube = anomalycube(magBP,cmBP,time,'event_dates.mat');

% apply land-sea-mask just in case
magBP_cube = landmaskcube(magBP_cube,lat,lon);

% average wind stress magnitude anomaly and 95% confidence interval mask
alpha = 0.05;
N = size(magBP_cube,3);
[magBP_avg_stat_sig,magBP_avg] = avg_anomaly_significance_mask(magBP_cube,N,alpha);

% mapping the average of the wind stress magnitude anomaly
titlestr = " ";
clabel = "$|\vec{\tau}|'\textsf{ [N m}^{-2}\textsf{]}$";
cmap = 'balance';
[~,~,~] = avgmap(magBP_avg_stat_sig,lat,lon,titlestr,clabel,cmap);
c = findall(gcf,'type','ColorBar');
c.Label.FontSize = 8;
plotm(y1,x1,'r','LineWidth',0.7)
plotm(y2,x2,'r','LineWidth',0.7)
plotm(q1,p1,'Color',[0.4940 0.1840 0.5560],'LineWidth',0.7)
c.FontSize = 8;
set(gca,'FontSize',8)
setm(gca,'FontSize',8)
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');
PL1.Children.ZData = 4;
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');
PL2.Children.ZData = 5;
t1 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8);

set(gcf,'Units','centimeters','Position',[0 0 9.5 11],'PaperUnits','centimeters','PaperPosition',[0 0 9.5 11])
% exportgraphics(gcf,'avg_windstr_magBP_stat_sig_v8.pdf','Resolution',1600)

% saveas(gcf,'avg_windstr_magBP_stat_sig_v6.png')

%% 14) Individual event maps with bandpass filtered wind stress magnitude anomaly
% grouptitle = "Daily 10-day to 6-month bandpass filtered wind stress magnitude anomaly";
% clabel = "$$ |\vec{\tau}|$$ [N/m$^{-2}$]";
% eventmap(magBP_cube,'event_dates.mat',lat,lon,grouptitle,clabel)
% 
% saveas(figure(2),'ev_windstr_magBP_1_v1.png')
% saveas(figure(3),'ev_windstr_magBP_2_v1.png')
% saveas(figure(4),'ev_windstr_magBP_3_v1.png')
% 
% %% 15) LP composite magnitude anomaly map
% magLP_cube = anomalycube(magLP,cmLP,time,'event_dates.mat');
% 
% % apply land-sea-mask just in case
% magLP_cube = landmaskcube(magLP_cube,lat,lon);
% 
% % average wind stress magnitude anomaly and 95% confidence interval mask
% alpha = 0.05;
% N = size(magLP_cube,3);
% [magLP_avg_stat_sig,magLP_avg] = avgstatsigmask(magLP_cube,N,alpha);
% 
% % mapping the average wind stress magnitude anomaly
% titlestr = "Event average 10-day low-pass filtered wind stress magnitude anomaly";
% clabel = "$$ |\vec{\tau}|$$ [N/m$^{-2}$]";
% [~,~,~] = avgmap(magLP_avg_stat_sig,lat,lon,titlestr,clabel);
% 
% % saveas(gcf,'avg_windstr_magLP_stat_sig_v1.png')
% 
% %% 16) LP magnitude anomaly event maps
% grouptitle = "Daily 10-day low-pass filtered wind stress magnitude anomaly";
% clabel = "$$ |\vec{\tau}|$$ [N/m$^{-2}$]";
% eventmap(magLP_cube,'event_dates.mat',lat,lon,grouptitle,clabel)
% 
% % saveas(figure(6),'ev_windstr_magLP_1_v1.png')
% % saveas(figure(7),'ev_windstr_magLP_2_v1.png')
% % saveas(figure(8),'ev_windstr_magLP_3_v1.png')
% %% 17) Map of mean summer magnitude of wind stress, not band-pass filtered
% load windstrmag.mat
% 
% [~,M,~,~,~,~] = datevec(time);
% keep = M > 11 | M < 3;
% mag_summer = windclim(:,:,keep);
% 
% % apply land-sea-mask just in case
% mag_summer = landmaskcube(mag_summer,lat,lon);
% 
% mag_summean = mean(mag_summer,3);
% high = mag_summean >= 0.3;
% low = mag_summean <= 2*10^(-2);
% mag_summean(high) = 0.3;
% mag_summean(low) = 0;
% 
% % mapping the average of the wind stress magnitude anomaly
% titlestr = "Average summer wind stress magnitude";
% clabel = "$$ |\vec{\tau}|$$ [N/m$^{-2}$]";
% cmap = 'balance';
% [~,~,~] = avgmap(mag_summean,lat,lon,[],clabel,cmap);
% 
% % saveas(figure(5),'avg_windstr_windclim_summer_v3.png')
% 
% %% 18) overall mean wind stress magnitude on same time scales
% 
% mag_mean = mean(windstrmag,3,'omitnan');
% high = mag_mean >= 0.3;
% low = mag_mean <= 2*10^(-2);
% mag_mean(high) = 0.3;
% mag_mean(low) = 0;
% 
% clabel = "$$ |\vec{\tau}|$$ [N/m$^{-2}$]";
% cmap = 'balance';
% 
% [~,~,~] = avgmap(mag_mean,lat,lon,[],clabel,cmap);
% 
% % saveas(figure(7),'avg_windstr_windstrmag_annual_v1.png')

% %% wind stress curl
% % Kylene 2022-04-19: I think I ended up not using this for wind stress
% % curl/Ekman pumping velocity figures
% 
% % load QuikSCAT 11/1/1999-10/31/2009
% qscat = load('/home/eddie/data4/loneill/exp_output/mar22_exp7-output.mat');
% crlQraw = qscat.strcrl;
% crlQ = crlQraw(:,:,1:14:end);
% time = squeeze(qscat.time(1:14:end));
% 
% % daily average
% x = 1;
% for i = 1:length(time)
%     while x<length(time)
%     crlQ(:,:,i) = mean(crlQraw(:,:,x:x+13),3,'omitnan');
%     x = x+14;
%     end
% end
% 
% % daily anomalies and climatology
% climQ = zeros([size(crlQ,1,2) 366]);
% anomQ = zeros(size(crlQ));
% d = day(qscat.time,'dayofyear');
% %i = 0;
% %while i < length(crlQ,3)
% for j = 1:366
%     d0 = d==j;
%     climQ(:,:,j) = mean(crlQ(:,:,d0),3,'omitnan');
%     ind0 = find(d0);
%     for k = sum(d0)
%         ind = ind0(k);
%         %anom
%     end
% 
% end
% % load ASCAT-A 
% % KNMI 25 km from 2007 through 2009
% % 12 km coastal processing 9/1/2010 to 8/31/2021
% ascat12 = load('/home/eddie/data4/loneill/exp_output/mar22_exp4-output.mat');
% crlA12 = ascat12.strcrl;
%  
% for i = 1:14
%     strcrl = squeeze(crlA12(:,:,i));
% end
%     