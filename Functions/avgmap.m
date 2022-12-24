% avgmap.m

% Maps a provided 2D composite average in the statistically significant avg
% matrix input on a map projection using Mapping Toolbox commands. The 'ax'
% argin (input argument) for an existing axes object on which to plot is
% optional, and it is useful for plotting maps in a tiled layout. To omit
% the plot title, the colorbar label, or use the matlab default colormap
% Jet instead of a cmocean colormap, respectively titlestr, clabel, and
% cmomap may be replaced by an empty vector, [].

% Inputs:
% dat_avg_stat_sig: a provided 2D matrix, which must already have values to be
%       plotted or ignored manipulated outside of the function.
% lat: a vector containing all latitudes within the map bounds, inclusive
%       of endpoints.
% lon: a vector containing all longitudes within the map bounds, inclusive
%       of endpoints.
% titlestr: the title, written in LaTeX as a string object.
% clabel: the color bar label, written in LaTeX as a string object.
% cmomap: the name of the chosen colormap from the cmocean package
%       developed by Thyng et al. (2016) as a string object.
% ax: an optional axes object on which the map created in this function
%       will be plotted. 

% Outputs: 
% C_dat: contour matrix.
% h_dat: contour patches.
% c: the color bar object.

% Kylene Cooley
% 26 Aug 2021
% Edited 13 Nov 2022 by KMC

function [C_dat,h_dat,c] = avgmap(dat_avg_stat_sig,lat,lon,titlestr,clabel,cmomap,ax)
% Coastline data set and coordinate limits around Chile-Peru System:
latlim = [min(lat) max(lat)];                   % latitude limits in [southern, northern]
lonlim = [min(lon) max(lon)];                   % longitude limits [western, eastern], using min/max this way works as long as lon is in [0 360]
load coastlines                                 % dataset included in Mapping Toolbox
Lat = lat*ones(length(lon),1)';                 % 2D latitude matrix with one latitudinal value/row (outer product)
Lon = ones(length(lat),1)*lon';                 % 2D longitude matrix with one longitudinal value/column
clim = [min(dat_avg_stat_sig,[],'all',...       % limits for the colorbar from data matrix
    'omitnan') max(dat_avg_stat_sig,[],...
    'all','omitnan')]; 
lvls = [min(clim,[],'omitnan'):10^(-3):...      % vector of contour level values
    max(clim,[],'omitnan')]; 

if exist('ax','var')                            % when an axes object is provided...
    subplot(ax)                                 % create new subplot in those axes
else
    figure()                                    % otherwise create new figure
end

h = worldmap(latlim,lonlim);                    % add a map projection to axes with latitude and longitude limits set above
setm(h,'PLabelLocation',latlim,...              % set parallel and meridian label parameters 
    'MLabelLocation',lonlim-360,...
    'MLabelParallel','south') 
[C_dat,h_dat] = contourm(Lat,Lon,...            % plot contour of SST' and returns contour matrix and patches
    dat_avg_stat_sig,lvls,'Fill','on'); 
caxis(clim)                                     % set axis limits for the colorbar

% This next part uses the argin 'strcmp' to retrieve the specified
% color map from the cmocean package as described in Thyng et al.
% (2016). This package must be downloaded locally and unzipped in the
% MATLAB Path. It is available from the File Exchange or in the Add-On
% browser. By downloading and opening the package within MATLAB, I
% found the package in the directory called 'MATLAB Add-Ons'. In this
% function, only 'balance' will return a pivot of 0 since it was the only
% diverging colormap that I used on this project.
if strcmp(cmomap,'balance')                     % for the 'balance' color map:
    cmocean('balance',length(lvls),'pivot',0)   % call this colormap from cmocean with the same number of levels as the number of contours and set the center of the diverging colormap to equal zero
else                                            % for other color maps:
    cmocean(cmomap,length(lvls))                % call this colormap from cmocean with the same number of levels as the number of contours
end

c = colorbar();                                 % retrieve colorbar object from current axes
fillm(coastlat,coastlon,[0.7 0.7 0.7])          % add gray landmass
c.Label.Interpreter = 'latex';                  % set colorbar label interpreter to LaTeX
c.Label.String = clabel;                        % add clabel as colorbar label after changing interpreter to avoid warning
c.Label.FontSize = 8;                           % set colorbar label font size
title(titlestr,'Interpreter','tex','FontSize',8)% add title to axes with tex interpreter and same font size
end
