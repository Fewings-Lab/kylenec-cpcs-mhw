% avg_anomaly_significance_mask.m

% Applies a mask to the average of provided anomalies assuming that the
% data follows a Gaussian distribution and averages that fall within a
% confidence interval about zero are not statistically significant from
% zero. Therefore, those averages are not significantly different from the
% climatological annual cycle. For N*=N-1 degrees of freedom the function
% evaluates the (100-alpha/2)% confidence interval on the average over N
% observations in the 3D matrix dat_cube. The function then uses the
% confidence interval to determine where the average is outside of the
% confidence interval about zero.

% Inputs:
% dat_cube: 3D matrix of anomalies with size [K, M, N]
% N: scalar number of observations, size(dat_cube,3)
% alpha: scalar confidence level

% Outputs:
% dat_avg_stat_sig: 2D matrix of average anomalies with significance mask
%   applied, size [K, M]
% dat_avg: 2D matrix of average anomalies without mask

% Kylene Cooley
% 26 Aug 2021
% Edited 13 Nov 2022 by KMC

function [dat_avg_stat_sig,dat_avg] = avg_anomaly_significance_mask(dat_cube,N,alpha)
dat_avg = mean(dat_cube,3,'omitnan');

sigma = std(dat_cube,1,3,'omitnan');

p = 1-(alpha/2);
q_t = tinv(p,N-1);

delta_mu = (q_t/sqrt(N)).*sigma;

% mask of statistically significant means
sigmask = xor(dat_avg>delta_mu,dat_avg<-delta_mu);

dat_avg_stat_sig = dat_avg;
dat_avg_stat_sig(~sigmask) = NaN;
end