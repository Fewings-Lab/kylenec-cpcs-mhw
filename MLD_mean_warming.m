% MLD_mean_warming.m 
% Kylene Cooley
% 10 Feb 2022
% Use outputs of MLD linear regression (b1) and the summer mean MLD from
% Holte et al. 2017 Argo MLD climatology with the mean anomalous warming
% outline to find mean MLD of each map within the area of strong warming.
% Includes a confidence interval on each mean and the largest whole number
% factor difference between each mean (not propagating error bc this is
% just for a talk).
% Edited 2 Mar 2022 to weight summer mean by area of square as a function
% of latitude
% Edited 25 Jan 2023 to make data paths more general using data sub-folders

currentdirectory = pwd;

% Staring with the linear regression output to test our method. Holte et
% al. summer mean MLD is used towards end of code.
mldLinearRegressionFile = fullfile(currentdirectory,'data','processed','mld_linreg_out.mat');
load(mldLinearRegressionFile)
% contains the output terms of the MLD linear regression model assuming no
% residual in the simplified heat budget, skill of the model, and lat, lon
% pairs of each point: b0, b1, lat, lon, S, Scr
MeanAreaWarmingFile = fullfile(currentdirectory,'data','processed','0_05_level_lines.mat');
load(MeanAreaWarmingFile)
% contains contour lines for mean area of anomalous warming
% top line x2,y2 and bottom line x1,y1

% Mean anomalous warming area lines have multiple points per longitude 
%   - size(x1) = 1x127 and size(x2) = 1x107 doubles
%   - x's are longitude, y's are latitude
% The problem is that these longitude's don't line up exactly with the lon
% 81x1 vector.
% does x1 == x2? no, x1 and x2 have different lengths
% could interpolate or project lines onto lat,lon regular grid
% then these points on grid = 1 and elsewhere is zero
latQ1 = interp1(x1,y1,lon); % for the line x1,y1 we are querying what latitude is nearest y1 if we interpolate at the coordinates in the vector "lon".
% Let's do same for x2,y2 and map these two pairs of lines on same axes
latQ2 = interp1(x2,y2,lon); % same as above for line x2,y2 at query coordinates in lon

%% map of x1,y1 & x2,y2 & lon,latQ1 & lon,latQ2 on same map axes
% basic info for maps
latlim = [min(lat) max(lat)];
lonlim = [min(lon) max(lon)]; 
load coastlines

% open a fresh map, add some labels and coastlines
figure()
g = worldmap(latlim,lonlim);
setm(g,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')
fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass

% Our original contour in red
h(1) = plotm(y1,x1,'r','LineWidth',2,'DisplayName','Original');
plotm(y2,x2,'r','LineWidth',2)

% Interpolated contour in blue
h(2) = plotm(latQ1,lon,'b','LineWidth',2,'DisplayName','Interpolation');
plotm(latQ2,lon,'b','LineWidth',2)

% add a legend
legend(h)

set(gcf,'PaperPosition',[0 0 5 6])
% saveas(gcf,'0_05_level_interp.png')

%% Using interpolated latQ1&2, mask area between the contour lines.
% We have a top and bottom latitude for every longitude on our coordinate
% grid. Based on the original numbering, top lat is latQ2 and bottom is latQ1.

% start mask w/ array of ones
areamask = ones(size(b1));

% next part requires knowing how lines 1 and 2 look, line 1 ends before
% line 2 in longitude
% find index where lines y1 and y2 end
last1 = find(isnan(latQ1),1,'first')-1; 
last2 = find(isnan(latQ2),1,'first')-1;

% use interpolated lines to make latitudes outside the lines = 0
for i= 1:length(lon)
    if i>last1 % solution after longitude where line 1 ends
        bottom = latQ1(last1); % use the last latitude in the line
    else
        bottom = latQ1(i);
    end
    if i>last2 % after last longitude in line 2
        top = latQ2(last2); % same solution applied here
    else
        top = latQ2(i);
    end
    latmask = lat>=bottom & lat<=top; % save latitudes between top and bottom lines
    areamask(:,i) = latmask; % save to area mask
end

% apply land mask after latmask & use land mask to make land = 0
landmask = read_era5_new('other','daily','01/01/1979','12/31/2020',{'land_sea_mask'},lonlim,latlim);
boo = find(round(landmask.land_sea_mask));
areamask(boo) = 0;

%% map to check area
% open a fresh map, add some labels and coastlines
figure()
g = worldmap(latlim,lonlim);
setm(g,'PLabelLocation',latlim,'MLabelLocation',lonlim-360,'MLabelParallel','south')

surfacem(lat,lon,areamask)

% Our original contour in red
h(1) = plotm(y1,x1,'r','LineWidth',2,'DisplayName','Original');
plotm(y2,x2,'r','LineWidth',2)

% Interpolated contour in blue
h(2) = plotm(latQ1,lon,'b','LineWidth',2,'DisplayName','Interpolation');
plotm(latQ2,lon,'b','LineWidth',2)

fillm(coastlat,coastlon,[0.7 0.7 0.7]) % Adds gray landmass
set(gcf,'PaperPosition',[0 0 5 6])
% saveas(gcf,'0_05_level_mask.png')

%% Save mask variable
save('0_05_level_mask.mat','areamask','lat','lon','latQ2','latQ1','latlim','lonlim')

%% Load linear regression MLD output and find mean in mean peak warming area
% load file
load mld_linreg_out.mat
load 0_05_level_mask.mat % uncomment if starting here

b1warm = b1.*areamask; % apply area mask to the linear regression output
% need anything = 0 to be nan
zip = b1warm == 0;
b1warm(zip)=nan; % linear regression coefficients only within anomalous warming area
% b1mean = mean(b1warm,'all','omitnan'); % find mean

% weighted mean only within warming area
wi = cosd(lat).*isfinite(b1warm);
wiXi = wi.*b1warm;
W = sum(wiXi,'all','omitnan')/sum(wi,'all','omitnan');

% 95% confidence intervals on the mean (BUT NOT ON A WEIGHTED MEAN, CHECK
% THIS)
% sample variance of the parent population
w = sum(wi,'all','omitnan');
w2 = sum(wi.^2,'all','omitnan');
s2 = (w/(w^2-w2)) * sum(wi.*((b1warm-W).^2),'all','omitnan');

alpha = 0.05;

N_toomany = sum(b1warm>0,'all','omitnan');
N_weights = (w^2)/w2; % effectively the same as the number of points above but this is the N_eff for a weighted mean
N = N_toomany/16; %This is the greatest effective sample size possible since only the events are significant
sigma = sqrt(s2);

% sigma = std(b1warm,1,'all','omitnan');

p = 1-(alpha/2);
q_t = tinv(p,N-1);

delta_mu = (q_t/sqrt(N)).*sigma;

%% Load summer mean Holte et al. MLD output and find mean in mean peak warming area
% load file
load mld_mean_eval.mat
% this is tricky bc summer_mean is 36x21 double
% we will have to use a subset of the area mask
load 0_05_level_mask.mat % uncomment if starting here

summerwarm = summer_mean.*areamask(1:4:end,1:4:end); % apply area mask to the linear regression output
% need anything = 0 to be nan
zip = summerwarm == 0;
summerwarm(zip)=nan;

% Find a mean summer MLD and weight each grid square by area as a function of latitude
% After simplifying W = sum(weights*#)/sum(weights), weights are simply
% cosd(lat).
Llat = lat(1:4:end);

% weighted mean only within warming area
wi = cosd(Llat).*isfinite(summerwarm);
wiXi = wi.*summerwarm;
W = sum(wiXi,'all','omitnan')/sum(wi,'all','omitnan');

% sample variance of the parent population
w = sum(wi,'all','omitnan');
w2 = sum(wi.^2,'all','omitnan');
s2 = (w/(w^2-w2)) * sum(wi.*((summerwarm-W).^2),'all','omitnan');


% 95% confidence intervals on the mean
alpha = 0.05;

N_toomany = sum(summerwarm>0,'all','omitnan');
N = (w^2)/w2; % effectively the same as the number of points above but this is the N_eff for a weighted mean
sigma = sqrt(s2);

p = 1-(alpha/2);
q_t = tinv(p,N-1);

delta_mu = (q_t/sqrt(N)).*sigma;

%% Factor difference between means
fact = floor(summermeanmean/b1mean); % by using floor, I can say more than whatever factor this is difference