%%% EXPERIMENT EXPLANATION:
%%% By using a simple grid of a WLAN sharing 4 channels, we compute the
%%% performance achieved by applying e-greedy.

clc
clear all

disp('-----------------------')
disp('MMRS (Toy scenario)')
disp('-----------------------')

% Generate constants from 'constants.m'
constants

initialEpsilon = 1;

% Setup the scenario: generate WLANs and initialize states and actions
wlans = generate_network_3D(nWlans, 'line', 2, drawMap); % SAFE CONFIGURATION

[ tx_power_SR_MMRS, tx_power_RD_MMRS, instantaoiAfterAction_A_MMRS, instantaoiAfterAction_B_MMRS, time1, time2] = benchmarkAoI_MMRS(wlans);

save('MMRS_RS_(MultiFrequency)_0327.mat');
% [aoi_benchmark_MS] = benchmarkAoI_MS(wlans);
% 
% Plot the results
% if plotResults
%     display_results_individual_performance(aoi_evolution_sequential_eg, aoi_benchmark_MMRS, aoi_benchmark_MS,...
%         tx_power_SR_sequential_eg, tx_power_RD_sequential_eg, times_arm_has_been_played_sequential_eg);
% 
% end
