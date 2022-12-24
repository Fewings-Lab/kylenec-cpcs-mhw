% sensitivity_test.m

% For dependence on Q'net and R on MLD h. Is the h needed to make R
% disappear reasonable given what we know about the coastal dynamics?

% open mat files containing Q'net and dSST'/dt, and event dates if needed
% for titles
Qnet = load('qnet_ev_data.mat');
dSST = load('dSSTdt_cube.mat');
load event_dates.mat

% cp = 3850; % J/kg/degC
% rhow = 1025; % kg/m^3
h0 = 25; % m, the estimate of MLD used in temperature change from Q'net

% If R = 0, rearrange eq. # to find:
% h = Qnet'/(rhow*cp)*(dSST'dt)^(-1)

dT_Qn = Qnet.dat_cube_dT .* h0; % (Qnet./(rhow*cp))
dT_Qn_avg = mean(dT_Qn,3,'omitnan');
dSST_avg = mean(dSST.dat_cube,3,'omitnan');
h_avg = dT_Qn_avg.*((dSST_avg).^(-1));

% Evaluate the average of the MLD data cube and 95% CI at each location
% h_avg = mean(h,3,'omitnan');
% alpha = 0.05;
% N = size(h,3);
% sigma = std(h,1,3,'omitnan');
% p = 1-(alpha/2);
% q_t = tinv(p,N-1);
% delta_mu = (q_t/sqrt(N)).*sigma;

% mask to change values to 0 if 95% confidence interval has
% negative values
% sigmask = h_avg+delta_mu<0;

% h_avg_stat_sig = h_avg;
% h_avg_stat_sig(sigmask) = NaN;
% 
% zedmask = and(h_avg+delta_mu>0,h_avg-delta_mu<0);
% h_avg_stat_sig(zedmask) = 0;

% apply land mask
latlim = [min(Qnet.info.sshf.lat) max(Qnet.info.sshf.lat)];
lonlim = [min(Qnet.info.sshf.lon) max(Qnet.info.sshf.lon)]; % works for min/max if lon is in [0 360]
em = read_era5_new('other','daily','01/01/1979','12/31/2020',{'land_sea_mask'},lonlim,latlim);
boo = find(round(em.land_sea_mask));
h_avg(boo) = NaN;


% Coastline data set and coordinate limits around Chile-Peru System:

load coastlines
Lat = Qnet.info.sshf.lat*ones(length(Qnet.info.sshf.lon),1)';
Lon = ones(length(Qnet.info.sshf.lat),1)*Qnet.info.sshf.lon';
clim = [0 50]; % limits for the colorbar

% Map of MLD needed to make R disappear
figure()
g = worldmap(latlim,lonlim);
setm(g,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,h_avg,[0:0.5:50],'Fill','on'); % Contour of SST'
caxis(clim)
cmocean('thermal',100)
c = colorbar();
c.Label.Interpreter = 'latex';
c.Label.String = "h [$m$]";
c.Label.FontSize = 18;
sgtitle("Mean mixed layer depth if $ R = 0 $",'Interpreter','latex','FontSize',20)

%% MLD statistics
Mx = max(h_avg,[],'all');
Mn = min(h_avg,[],'all');

mu = mean(h_avg,'all','omitnan');
alpha = 0.05;
N = sum(isfinite(h_avg),'all');
sigma = std(h_avg,1,'all','omitnan');
p = 1-(alpha/2);
q_t = tinv(p,N-1);
delta_mu = (q_t/sqrt(N)).*sigma;

