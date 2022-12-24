% nearshoreboxmean.m
% A function to slice the data in the box where we found SST' events, and
% then take the mean of it to obtain a time series. Currently will set up
% to have 3 data cube inputs and 3 time series outputs. want to make this 
% optional if possible.
% Kylene Cooley
% 23 Aug 2021

function [ts1, ts2, ts3] = nearshoreboxmean(lat,lon,cube1,cube2,cube3)
    % event box coordinates
    latlim = [-36.5 -35.5];
    lonlim = 360 + [-74.25 -73.25];

    blat = find(sum(lat == latlim,2));
    blon = find(sum(lon == lonlim,2));

    foo1_1 = cube1(blat(1):blat(2),blon(1):blon(2),:);
    ts1 = squeeze(mean(foo1_1,[1 2],'omitnan'));
    
    foo1_2 = cube2(blat(1):blat(2),blon(1):blon(2),:);
    ts2 = squeeze(mean(foo1_2,[1 2],'omitnan'));

    foo1_3 = cube3(blat(1):blat(2),blon(1):blon(2),:);
    ts3 = squeeze(mean(foo1_3,[1 2],'omitnan'));
end