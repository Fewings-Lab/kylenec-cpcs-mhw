% s07_2_MLD_regression.m

% Working through process of cross-correlation map between dSST' and
% Q'_net/rhow_w*c_p, to linear regression between these to find the slope
% of the line of best fit, h

% Edited 2 Mar 2022 to make a figure with linear regression and summer
% climatological MLD as subplots

% Edited 13 Nov 2022 with comments by KMC

% open mat files containing Q'net and dSST'/dt
Qnet = load('qnet_ev_data.mat');
dSST = load('dSSTdt_cube.mat');

% cp = 3850; % J/kg/degC
% rhow = 1025; % kg/m^3
h0 = 25; % m, the estimate of MLD used in temperature change from Q'net
% If R = 0, rearrange eq. # to find:
% h = Qnet'/(rhow*cp)*(dSST'dt)^(-1)

dT_Qn = Qnet.dat_cube_dT .* h0; % (Qnet./(rhow*cp))
dSSTdt = dSST.dat_cube;

%% Check if dSST' and Q'_net/rho_w*c_p is appears correlated-ish

% Picked location based off h map from sensitivity test
Qpt = dT_Qn(101,38,:);
Tpt = dSSTdt(101,38,:);

% Scatter plot with y1 v. y2
figure()
scatter(Qpt,Tpt,4,'filled')
title("Scatter plot of $\frac{\partial SST'}{\partial t}$ against $ \frac{Q'_{net}}{\rho_w c_p}$",'Interpreter','latex')
xlabel("$ \frac{Q'_{net}}{\rho_w c_p}$",'Interpreter','latex')
ylabel("$\frac{\partial SST'}{\partial t}$",'Interpreter','latex')

%% Correlation coefficient map

rhoSSTqnet = zeros(size(dSSTdt,[1 2]));
alpha = 0.05; % 95% confidence level
Nstr = 37;
rho_cr = (finv(1-alpha,1,Nstr-2)/(Nstr-2+finv(1-alpha,1,Nstr-2)))^(1/2); % using Fisher-F dist.

% loop through lat and lon
for lat = 1:size(dSSTdt,1)
    for lon = 1:size(dSSTdt,2)
        xx = dT_Qn(lat,lon,:);
        yy = dSSTdt(lat,lon,:);
        
        Rxy= acov(xx,yy);
        
        rhoxy = acor(Rxy,xx,yy);
        
        rhoSSTqnet(lat,lon)= rhoxy;
    end
end

boo = rhoSSTqnet < rho_cr;
rhoSSTqnet(boo) = NaN;

% apply land mask
lat = Qnet.info.sshf.lat;
lon = Qnet.info.sshf.lon;
latlim = [min(Qnet.info.sshf.lat) max(Qnet.info.sshf.lat)];
lonlim = [min(Qnet.info.sshf.lon) max(Qnet.info.sshf.lon)]; % works for min/max if lon is in [0 360]
em = read_era5_new('other','daily','01/01/1979','12/31/2020',{'land_sea_mask'},lonlim,latlim);
boo = find(round(em.land_sea_mask));
rhoSSTqnet(boo) = 0;

% % Coastline data set and coordinate limits around Chile-Peru System:
% load coastlines
% Lat = Qnet.info.sshf.lat*ones(length(Qnet.info.sshf.lon),1)';
% Lon = ones(length(Qnet.info.sshf.lat),1)*Qnet.info.sshf.lon';
% clim = [min(rhoSSTqnet,[],'all','omitnan') 1]; % limits for the colorbar
% lvls = [min(clim):0.005:max(clim)];
% 
% % map correlation coefficient rhoSSTqnet
% figure(8)
% g = worldmap(latlim,lonlim);
% setm(g,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
% plotm(coastlat,coastlon) % Adds coastlines
% [C,~] = contourm(Lat,Lon,rhoSSTqnet,lvls,'Fill','on'); % Contour of SST'
% caxis(clim)
% cmocean('balance',length(lvls),'pivot',0)
% % cmocean('balance',length(lvls))
% c = colorbar();
% set(c,'FontSize',14)
% c.Label.Interpreter = 'latex';
% c.Label.String = "$\rho_{\frac{Q'_{net}}{\rho_w c_p},\frac{\partial SST'}{\partial t}}$";
% setm(g,'FontSize',14)
% set(c.Label,'FontSize',24)
% % sgtitle("Correlation coefficient $\rho$ between $\frac{Q'_{net}}{\rho_w c_p}$ and $\frac{\partial SST'}{\partial t}$",'Interpreter','latex','FontSize',24)
cmap = 'balance';
clabel = "$\hat{\rho}$";
[~,~,c] = avgmap(rhoSSTqnet,lat,lon,[],clabel,cmap);

c.FontSize = 8;
c.Label.FontSize = 12;
set(gcf,'Units','centimeters','Position',[0 0 9.5 11],'PaperUnits','centimeters','PaperPosition',[0 0 9.5 11])
set(gca,'FontSize',8)
setm(gca,'FontSize',8)
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');
PL1.Children.ZData = 4;
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');
PL2.Children.ZData = 5;
t1 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8);

% exportgraphics(gcf,'corr_coeff_Qnet_dSST_v8.pdf','Resolution',1600)


%% Linear regression, where MLD, h, is b1 coefficient
% first at one point (same as test above):
% 37 years from 1980 to 2019
y = squeeze(Qpt); % [y] = degC*m/day, Q'net/rho_w*c_p
[t,I] = sort(squeeze(Tpt)); % [t] = degC/day, dSST'/dt
y = y(I);

% Two-parameter, univariate regression y1 = b0+b1*t
% Using linear regression methods from OC 682 where b1 is MLD
Ryt = acov(y,t);
vart = var(t,1,'omitnan');
b1 = Ryt/vart;
muy = mean(y,'omitnan');
mut = mean(t,'omitnan');
b0 = muy-(Ryt*mut/vart);
y1 = b0 + b1.*t;
S = (acor(Ryt,y,t))^2;
Nst = 37;

% Confidence intervals on regression coefficients
Y = y;
x1 = t;
NotNaN = ~isnan(Y.*x1);
N = sum(NotNaN);
Y = Y(NotNaN);
x1 = x1(NotNaN);
Z = [sum(Y,'omitnan')/N; sum(Y.*x1,'omitnan')/N];
D = [1, sum(x1,'omitnan')/N; sum(x1,'omitnan')/N,sum(x1.^2,'omitnan')/N];
B = D\Z;
%y1 = B(1) + B(2).*x1;

RYy1 = sum((Y-mean(Y,'omitnan')).*(y1-mean(y1,'omitnan')),'omitnan')/N;
varY = sum((Y-sum(Y,'omitnan')/N).^2,'omitnan')/N;
vary1 = sum((y1-sum(y1,'omitnan')/N).^2,'omitnan')/N;
S = (RYy1^2)/(varY*vary1);
M = 1;
alpha = 0.05;
% Confidence intervals on regression coefficients 
RBB = (((1-S)*varY)/(Nst-M-1)).*inv(D); % Covariance matrix of B with B
BCon = [sqrt(RBB(1,1));sqrt(RBB(2,2))].*tinv(1-(alpha/2),Nst-M-1);

% Critical skill based on a 95% significance test
Scr = (M*finv(1-alpha,M,Nst-M-1))/(Nst-M-1+M*finv(1-alpha,M,Nst-M-1));

% Plot y v. t and regression model, y1 v. t
% Note values of b0, b1 and skill S
figure(3)
plot(t,y,'.r',t,y1,'-k')
title('Mixed layer depth linear regression at point 101,38')
xlabel("$\frac{\partial SST'}{\partial t}$ [$ ^\circ $C/day]",'Interpreter','latex')
ylabel("$\frac{Q'_{net}}{\rho_w c_p}$ [$ ^\circ $C $ \times $ m/day]",'Interpreter','latex')
text(-0.1,1.0,['$\hat{S}_{crit}  =$ ',num2str(Scr,3)],'interpreter','latex')
text(-0.1,1.2,['$\hat{S}  =$ ',num2str(S,3)],'interpreter','latex')
text(-0.1,1.4,['h = \beta_1 = ',num2str(b1,3),'\pm',num2str(BCon(2),3)])
text(-0.1,1.6,['\beta_0 = ',num2str(b0,3),'\pm',num2str(BCon(1),3)])


%% then using all points
% 37 years from 1980 to 2019
y = reshape(dT_Qn,[numel(dT_Qn) 1]); % [y] = degC*m/day, Q'net/rho_w*c_p
[t,I] = sort(reshape(dSSTdt,[numel(dSSTdt) 1])); % [t] = degC/day, dSST'/dt
y = y(I);

% Two-parameter, univariate regression y1 = b0+b1*t
% Using linear regression methods from OC 682 where b1 is MLD
Ryt = acov(y,t);
vart = var(t,1,'omitnan');
b1 = Ryt/vart;
muy = mean(y,'omitnan');
mut = mean(t,'omitnan');
b0 = muy-(Ryt*mut/vart);
y1 = b0 + b1.*t;
S = (acor(Ryt,y,t))^2;
Nst = 37;

% Confidence intervals on regression coefficients
Y = y;
x1 = t;
NotNaN = ~isnan(Y.*x1);
N = sum(NotNaN);
Y = Y(NotNaN);
x1 = x1(NotNaN);
Z = [sum(Y,'omitnan')/N; sum(Y.*x1,'omitnan')/N];
D = [1, sum(x1,'omitnan')/N; sum(x1,'omitnan')/N,sum(x1.^2,'omitnan')/N];
B = D\Z;
%y1 = B(1) + B(2).*x1;
y1_sm = y1(NotNaN);

RYy1 = sum((Y-mean(Y,'omitnan')).*(y1_sm-mean(y1_sm,'omitnan')),'omitnan')/N;
varY = sum((Y-sum(Y,'omitnan')/N).^2,'omitnan')/N;
vary1 = sum((y1_sm-sum(y1_sm,'omitnan')/N).^2,'omitnan')/N;
S = (RYy1^2)/(varY*vary1);
M = 1;
alpha = 0.05;
% Confidence intervals on regression coefficients 
RBB = (((1-S)*varY)/(Nst-M-1)).*inv(D); % Covariance matrix of B with B
BCon = [sqrt(RBB(1,1));sqrt(RBB(2,2))].*tinv(1-(alpha/2),Nst-M-1);

% Critical skill based on a 95% significance test
Scr = (M*finv(1-alpha,M,Nst-M-1))/(Nst-M-1+M*finv(1-alpha,M,Nst-M-1));

% Plot y v. t and regression model, y1 v. t
% Note values of b0, b1 and skill S
figure()
plot(t,y,'.r',t,y1,'-k')
title('Mixed layer depth linear regression for all locations')
xlabel("$\frac{\partial SST'}{\partial t}$ [$ ^\circ $C/day]",'Interpreter','latex')
ylabel("$\frac{Q'_{net}}{\rho_w c_p}$ [$ ^\circ $C $ \times $ m/day]",'Interpreter','latex')
text(-0.6,2.1,['$\hat{S}_{crit}  =$ ',num2str(Scr,3)],'interpreter','latex')
text(-0.6,2.8,['$\hat{S}  =$ ',num2str(S,3)],'interpreter','latex')
text(-0.6,3.5,['h = \beta_1 = ',num2str(b1,3),'\pm',num2str(BCon(2),3)])
text(-0.6,4.2,['\beta_0 = ',num2str(b0,3),'\pm',num2str(BCon(1),3)])

%% then map at each location
b0 = zeros(size(dSSTdt,[1 2]));
b1 = zeros(size(dSSTdt,[1 2]));
S = zeros(size(dSSTdt,[1 2]));
Scr = zeros(size(dSSTdt,[1 2]));
db0 = zeros(size(dSSTdt,[1 2]));
db1 = zeros(size(dSSTdt,[1 2]));


% loop through lat and lon
for i = 1:size(dSSTdt,1)
    for j = 1:size(dSSTdt,2)
        y = squeeze(dT_Qn(i,j,:)); % [y] = degC*m/day, Q'net/rho_w*c_p
        t = squeeze(dSSTdt(i,j,:)); % [t] = degC/day, dSST'/dt

        % Two-parameter, univariate regression y1 = b0+b1*t
        % Using linear regression methods from OC 682 where b1 is MLD
        Ryt = acov(y,t);
        vart = var(t,1,'omitnan');
        b1(i,j) = Ryt/vart;
        muy = mean(y,'omitnan');
        mut = mean(t,'omitnan');
        b0(i,j) = muy-(Ryt*mut/vart);
        y1 = b0(i,j) + b1(i,j).*t;
        S(i,j) = (acor(Ryt,y,t))^2;
        Nst = 37;
        
        % Confidence intervals on regression coefficients
        Y = y;
        x1 = t;
        NotNaN = ~isnan(Y.*x1);
        N = sum(NotNaN);
        Y = Y(NotNaN);
        x1 = x1(NotNaN);
        Z = [sum(Y,'omitnan')/N; sum(Y.*x1,'omitnan')/N];
        D = [1, sum(x1,'omitnan')/N; sum(x1,'omitnan')/N,sum(x1.^2,'omitnan')/N];
        B = D\Z;
        %y1 = B(1) + B(2).*x1;
        y1_sm = y1(NotNaN);

        RYy1 = sum((Y-mean(Y,'omitnan')).*(y1_sm-mean(y1_sm,'omitnan')),'omitnan')/N;
        varY = sum((Y-sum(Y,'omitnan')/N).^2,'omitnan')/N;
        vary1 = sum((y1_sm-sum(y1_sm,'omitnan')/N).^2,'omitnan')/N;
        S(i,j) = (RYy1^2)/(varY*vary1);
        M = 1;
        alpha = 0.05;
        % Confidence intervals on regression coefficients 
        RBB = (((1-S(i,j))*varY)/(Nst-M-1)).*inv(D); % Covariance matrix of B with B
        BCon = [sqrt(RBB(1,1));sqrt(RBB(2,2))].*tinv(1-(alpha/2),Nst-M-1);
        db0(i,j) = BCon(1);
        db1(i,j) = BCon(2);
        
        % Critical skill based on a 95% significance test
        Scr(i,j) = (M*finv(1-alpha,M,Nst-M-1))/(Nst-M-1+M*finv(1-alpha,M,Nst-M-1));

    end
end

% critical skill mask
boo = S < Scr;
b0(boo) = NaN;
b1(boo) = NaN;
db0(boo) = NaN;
db1(boo) = NaN;

% apply land mask
latlim = [min(Qnet.info.sshf.lat) max(Qnet.info.sshf.lat)];
lonlim = [min(Qnet.info.sshf.lon) max(Qnet.info.sshf.lon)]; % works for min/max if lon is in [0 360]
em = read_era5_new('other','daily','01/01/1979','12/31/2020',{'land_sea_mask'},lonlim,latlim);
boo = find(round(em.land_sea_mask));
b0(boo) = 0;
b1(boo) = 0;
db0(boo) = 0;
db1(boo) = 0;

%% Map MLD regression and summer mean MLD with large cbar range
% Coastline data set and coordinate limits around Chile-Peru System:
lat = Qnet.info.sshf.lat;
lon = Qnet.info.sshf.lon;

load coastlines
Lat = lat;
Lon = lon;

load('mld_mean_eval.mat','summer_mean')
clim = [0 max(summer_mean,[],'all','omitnan')]; % limits for the colorbar
lvls = [min(clim):0.005:max(clim)];

% map correlation coefficient rhoSSTqnet
figure()
tiledlayout(1,2,'TileSpacing','tight','Padding','loose')
% what would 'TileSpacing' 'loose' and 'Padding' 'tight' look like?

nexttile
%subplot(1,2,1)
g = worldmap(latlim,lonlim);
setm(g,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
[C,~] = contourm(Lat,Lon,b1,lvls,'Fill','on'); % Contour of SST'
caxis(clim)
% cmocean('balance',length(lvls),'pivot',0)
cmocean('deep',length(lvls))
c = colorbar();
% set(c,'FontSize',14)
c.Label.Interpreter = 'latex';
c.Label.String = "$\textsf{MLD, }\hat{h}\textsf{ [m]}$";
% setm(g,'FontSize',14)
set(c.Label,'FontSize',8)
% sgtitle("Mixed layer depth by linear regression between $\frac{Q'_{net}}{\rho_w c_p}$ and $\frac{\partial SST'}{\partial t}$",'Interpreter','latex','FontSize',24)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
load('0_05_level_lines.mat')
plotm(y1,x1,'r','LineWidth',0.7)
plotm(y2,x2,'r','LineWidth',0.7)
c.FontSize = 8;
set(gca,'FontSize',8)
setm(gca,'FontSize',8)
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');
PL1.Children.ZData = 4;
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');
PL2.Children.ZData = 5;
t1 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8);

nexttile
%subplot(1,2,2)
g = worldmap(latlim,lonlim);
setm(g,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
plotm(coastlat,coastlon) % Adds coastlines
pcolorm(Lat(1:4:end),Lon(1:4:end),summer_mean); % Contour of SST'
caxis(clim)
% cmocean('balance',length(lvls),'pivot',0)
cmocean('deep',length(lvls))
c = colorbar();
% set(c,'FontSize',14)
c.Label.Interpreter = 'latex';
c.Label.String = "$\textsf{MLD, }h\textsf{ [m]}$";
% setm(g,'FontSize',14)
set(c.Label,'FontSize',8)
% sgtitle("Mixed layer depth by linear regression between $\frac{Q'_{net}}{\rho_w c_p}$ and $\frac{\partial SST'}{\partial t}$",'Interpreter','latex','FontSize',24)
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
load('0_05_level_lines.mat')
plotm(y1,x1,'r','LineWidth',0.7)
plotm(y2,x2,'r','LineWidth',0.7)
c.FontSize = 8;
set(gca,'FontSize',8)
setm(gca,'FontSize',8)
PL1 = scatterm(-37.3,-73.3,40,'ow','filled');
PL1.Children.ZData = 4;
PL2 = scatterm(-37.3,-73.3,17,'k','filled','MarkerEdgeColor','w');
PL2.Children.ZData = 5;
t2 = text(610000,-4250000,'PL','Color',[0 0 0],'HorizontalAlignment','left','FontSize',8);


set(gcf,'Units','centimeters','Position',[0 0 19 11],'PaperUnits','centimeters','PaperPosition',[0 0 19 11])
% exportgraphics(gcf,'compare_MLD_v3.pdf','Resolution',1600)

%% map MLD regression for small cbar range
lat = Qnet.info.sshf.lat;
lon = Qnet.info.sshf.lon;
clabel = "MLD, $\hat{h}$ [m]";
cmap = 'deep';
[~,~,~] = avgmap(b1,lat,lon,[],clabel,cmap);
set(gcf,'PaperPosition',[0 0 5 6])

% do this next line manually after getting size right:
% saveas(gcf,'linear_regr_MLD_v8.png')

%% Find scatter points within ID box
latlim_box = [-36.5 -35.5];
lonlim_box = 360 + [-74.25 -73.25];


        
        
%% Functions
function R = crcov(xx,yy)
    % Cross-covariance (in time) between two signals xx and yy
    S = (xx-mean(xx,'omitnan')).*(yy-mean(yy,'omitnan'));
    R = sum(S,'omitnan')/sum(~isnan(S)); % Normalized by N-1
end

function rho = crcor(Rxy,xx,yy)
    % Cross-correlation in time between xx and yy with cross-covariance Rxy
    rho = Rxy/(std(xx,1,'omitnan')*std(yy,1,'omitnan'));
end
