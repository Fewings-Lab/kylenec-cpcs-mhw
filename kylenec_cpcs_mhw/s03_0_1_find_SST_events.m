% find_SST_events.m
% Kylene Cooley
% 29 June 2021

% 17 June 2021: Change color of orange bars in histogram to red.
% -- KC

latlim = [-36.5 -35.5];
lonlim = 360 + [-74.25 -73.25];

% open daily SST anomaly data struct within 1deg box
dat = read_era5_new('sfc','daily','01/01/1979','12/31/2020',{'sst'},lonlim,latlim,'anomaly',1);
fields = fieldnames(dat);
data = eval(join(['dat.',fields{8}],'.'));

% filter to time scales of 10 days to 6 months with pl66

% low-pass filter with a 10-day pl66 filter
dt = dat.time(2)-dat.time(1); % sampling interval in matdate (number of days)
Tl = 10; % half-amplitude cutoff period in days to match matdate
dat_low = NaN(size(data)); % this will be the 10-day low-pass filtered data cube
for i=1:length(dat.lon) % looping through longitude slices since pl66 only handles 2D matrices
    dat_low(i,:,:) = (pl66tn(squeeze(data(i,:,:)),dt*24,Tl*24))'; % pl66 assumes time in hours
end

% high-pass filter with 6-month pl66 filter
Th = hours(years(0.5));
dat_lolo = NaN(size(dat_low));
for j = 1:length(dat.lon)
    dat_lolo(j,:,:) = (pl66tn(squeeze(dat_low(j,:,:)),dt*24,Th))';
end

dat_band = dat_low-dat_lolo;
dat_band(:,:,1:2*round(Th/(dt*24))) = NaN;
dat_band(:,:,end-2*round(Th/(dt*24)):end) = NaN;

% take spatial mean of this box
dat_avg = squeeze(mean(dat_band,[1 2],'omitnan'));

% find standard deviation
sigma = std(dat_avg,'omitnan');

% find SST' greater than 2*std and local maxima
strong = dat_avg > (2*sigma);
localmax = islocalmax(dat_avg);

% limit to dec through feb
dv = datevec(dat.time);
monthnum = dv(:,2);
summer = monthnum<3 | monthnum==12;

% find events in SST anomaly
ind = strong & summer & localmax;
t_ev = dat.time(ind);

%% Redo SST anomaly histogram from paper with  same filtering as time series
allev = sum(strong & localmax);
fall = monthnum>2 & monthnum<6;
fallev = sum(strong & localmax & fall);
winter = monthnum>5 & monthnum<9;
winterev = sum(strong & localmax & winter);
spring = monthnum>8 & monthnum<12;
springev = sum(strong & localmax & spring);

ref_t = datetime(1,1:12,1);
ref_m = month(ref_t,'shortname');
% t = datetime(dv);
% m = month(t,'shortname');
weak = dat_avg<(-2*sigma);

% need to pad blu below with nans to same length as orang

% orang = categorical(monthnum(strong),1:12,ref_m);
x = sum(strong) - sum(weak);
pad = NaN(x,1);
C = categorical([monthnum(strong),[monthnum(weak);pad]],1:12,ref_m); % order is original orange first column, blue in second column
% 
% figure()
% histogram(C,ref_m)

f = figure();
hist(C,ref_m)
h = findobj(gca,'Type','bar');
h(1).FaceColor = [0 0.4470 0.7410]; % blue bars, seems reversed but bin counts add up to # weak
h(2).FaceColor = [1 0 0]; % red bars
h(2).DisplayName = 'over 2\sigma';
h(1).DisplayName = 'below -2\sigma';
legend()
ylabel('Total days 1979-2020')
xlabel('Month')

f.Units = 'centimeters';
f.Position = [0 0 9.5 4];
f.PaperUnits = 'centimeters';
f.PaperPosition = [0 0 9.5 4];
set(gca,'FontSize',8)

% exportgraphics(f,'high-low-month-anom-v12.pdf')
%% open dSST'/dt to find the maxima in warming before events within 1deg box
dat_dt = read_era5_new('sfc','daily','01/01/1979','12/31/2020',{'dsstdt'},lonlim,latlim,'anomaly',1);
fields_dt = fieldnames(dat_dt);
data_dt = eval(join(['dat_dt.',fields_dt{6}],'.'));

% filter dSST'/dt
dat_low_dt = NaN(size(data_dt)); % this will be the 10-day low-pass filtered data cube
for i=1:length(dat_dt.lon) % looping through longitude slices since pl66 only handles 2D matrices
    dat_low_dt(:,i,:) = (pl66tn(squeeze(data_dt(:,i,:)),dt*24,Tl*24))'; % pl66 assumes time in hours
end

% high-pass filter with 6-month pl66 filter
dat_lolo_dt = NaN(size(dat_low_dt));
for j = 1:length(dat_dt.lon)
    dat_lolo_dt(:,j,:) = (pl66tn(squeeze(dat_low_dt(:,j,:)),dt*24,Th))';
end

dat_band_dt = dat_low_dt-dat_lolo_dt;
dat_band_dt(:,:,1:2*round(Th/(dt*24))) = NaN;
dat_band_dt(:,:,end-2*round(Th/(dt*24)):end) = NaN;

% take spatial mean of this box
dat_avg_dt = squeeze(mean(dat_band_dt,[1 2],'omitnan'));

% identify peaks in dSST'/dt
peaks = islocalmax(dat_avg_dt);

% make list of dates of peaks in dat_avg_dt that are before each SST' event
t_dt = NaN(length(t_ev),1);
for k = 1:length(t_ev)
   ktimes = dat_dt.time < t_ev(k);
   kind = find(ktimes & peaks,1,'last');
   if sum(kind)==0
       t_dt(k) = NaN;
   else
       t_dt(k) = dat_dt.time(kind);
   end
end


% save matdate times to .mat file
save('event_dates.mat','t_ev','t_dt')