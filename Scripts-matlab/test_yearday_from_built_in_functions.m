% test new simpler method for calculating yearday using only built-in
% Matlab functions

dn = datenum(2000,1,1:1/24:1000).'; % a few years of hourly time variable

datestr(dn(1:10)) % check contents of dn by converting to strings
% ans =
% 
%   10×20 char array
% 
%     '01-Jan-2000 00:00:00'
%     '01-Jan-2000 01:00:00'
%     '01-Jan-2000 02:00:00'
%     '01-Jan-2000 03:00:00'
%     '01-Jan-2000 04:00:00'
%     '01-Jan-2000 05:00:00'
%     '01-Jan-2000 06:00:00'
%     '01-Jan-2000 07:00:00'
%     '01-Jan-2000 08:00:00'
%     '01-Jan-2000 09:00:00'

datestr(dn(end-10:end)) % check contents of dn by converting to strings
% ans =
% 
%   11×20 char array
% 
%     '25-Sep-2002 14:00:00'
%     '25-Sep-2002 15:00:00'
%     '25-Sep-2002 16:00:00'
%     '25-Sep-2002 17:00:00'
%     '25-Sep-2002 18:00:00'
%     '25-Sep-2002 19:00:00'
%     '25-Sep-2002 20:00:00'
%     '25-Sep-2002 21:00:00'
%     '25-Sep-2002 22:00:00'
%     '25-Sep-2002 23:00:00'
%     '26-Sep-2002 00:00:00'

dv = datevec(dn); % convert to date vector [yy mm dd hh mm ss]
dv(1:10,:) % take a look
% ans =
% 
%         2000           1           1           0           0           0
%         2000           1           1           1           0           0
%         2000           1           1           2           0           0
%         2000           1           1           3           0           0
%         2000           1           1           4           0           0
%         2000           1           1           5           0           0
%         2000           1           1           6           0           0
%         2000           1           1           7           0           0
%         2000           1           1           8           0           0
%         2000           1           1           9           0           0

% calculate yearday by subtracting the datenum of Jan 1 of that year 
% (and add 1 so the first day of the year is yearday 1 not 0)
yd = dn - datenum(dv(:,1),1,1) + 1; 

% check whether this gives the same answer as mf_dn2yd
yd(end-10:end)
% ans =
% 
%   268.5833
%   268.6250
%   268.6667
%   268.7083
%   268.7500
%   268.7917
%   268.8333
%   268.8750
%   268.9167
%   268.9583
%   269.0000

mf_dn2yd(dn(end-10:end))
% ans =
% 
%   268.5833
%   268.6250
%   268.6667
%   268.7083
%   268.7500
%   268.7917
%   268.8333
%   268.8750
%   268.9167
%   268.9583
%   269.0000

% in summary:
dn = datenum(2000,1,1:1/24:1000).'; % a few years of hourly time variable
dv = datevec(dn); % convert to date vector [yy mm dd hh mm ss]
% calculate yearday by subtracting the datenum of Jan 1 of that year 
% (and add 1 so the first day of the year is yearday 1 not 0)
yd = dn - datenum(dv(:,1),1,1) + 1;

%%
figure(1)
clf
plot(mf_dn2yd(dn),yd,'r.')
refline(1,0)