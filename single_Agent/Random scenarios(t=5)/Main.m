%%% EXPERIMENT EXPLANATION:
%%% By using a simple grid of a WLAN sharing 4 channels, we compute the
%%% performance achieved by applying e-greedy.

clc
clear all

disp('-----------------------')
disp('e-greedy (Toy scenario)')
disp('-----------------------')

% Generate constants from 'constants.m'
constants
% Define EXP3 parameters
initialEta = 0.1;
gamma = 0;
initialEpsilon = 0.02;

% Setup the scenario: generate WLANs and initialize states and actions
wlans = generate_network_3D(nWlans, 'line', 2, drawMap); % SAFE CONFIGURATION

% Compute the AoI experienced per WLAN at each iteration
[aoi_evolution_concurrent_eg, times_arm_has_been_played_concurrent_eg, tx_power_SR_concurrent_eg, tx_power_RD_concurrent_eg, instantaoiAfterAction_A_concurrent_eg, instantaoiAfterAction_B_concurrent_eg] = ...
    concurrent_egreedy( wlans, initialEpsilon );
% [ aoi_evolution_concurrent_exp3, times_arm_has_been_played_concurrent_exp3, tx_power_SR_concurrent_exp3, tx_power_RD_concurrent_exp3, instantaoiAfterAction_A_concurrent_exp3, instantaoiAfterAction_B_concurrent_exp3] = ...
%     concurrent_exp3( wlans, gamma, initialEta );
% [ aoi_evolution_concurrent_ucb, times_arm_has_been_played_concurrent_ucb, tx_power_SR_concurrent_ucb, tx_power_RD_concurrent_ucb, instantaoiAfterAction_A_concurrent_ucb, instantaoiAfterAction_B_concurrent_ucb] = ...
%     concurrent_ucb( wlans );
% [ aoi_evolution_concurrent_ts, times_arm_has_been_played_concurrent_ts, tx_power_SR_concurrent_ts, tx_power_RD_concurrent_ts, instantaoiAfterAction_A_concurrent_ts, instantaoiAfterAction_B_concurrent_ts] = ...
%     concurrent_thompson_sampling( wlans );

% Save constants into current folder
save('single-Agent_RS_egreedy=0.02(t=5)10e5_(MultiFrequency)_0327.mat');

% Plot the results
% if plotResults
%     display_results_individual_performance(aoi_evolution_sequential_eg, aoi_benchmark_MMRS, aoi_benchmark_MS,...
%         tx_power_SR_sequential_eg, tx_power_RD_sequential_eg, times_arm_has_been_played_sequential_eg);
% 
% end
