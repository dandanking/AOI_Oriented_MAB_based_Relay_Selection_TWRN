nWlans = 4;                         % Number of WLANs in the map
totalIterations = 100000;            % Maximum convergence time (one period implies the participation of all WLANs)
minimumIterationToConsider = 5000;  % Iteration from which to consider the obtained results
totalScenarios = 100;               % Number of TOTAL repetitions to take the average

% Define the transitory and the permanent intervals of iterations
transitoryInterval = 1 : minimumIterationToConsider;
permanentInterval = minimumIterationToConsider + 1 : totalIterations;

NOISE_DBM = -10;                   % Floor NOISE_DBM in dBm

nChannels = 4;                      % Number of available channels (from 1 to n_channels)
channelActions = 1 : nChannels;     % Possible channels
ccaActions = [-82];                 % CCA levels (dBm)
txPowerActions = [-1 0 1 2];     % Transmit power levels (dBm)
benchmark_AOI = -5;

% Each state represents an [i,j,k] combination for indexes on "channels", "cca" and "tx_power"
possibleActions = 1:(size(channelActions, 2) * ...
    size(ccaActions, 2) * size(txPowerActions, 2));
K = size(possibleActions,2);   % Total number of actions
allCombs = allcomb(1:K, 1:K);

% Structured array with all the combinations (for computing the optimal)
possibleComb = allcomb(possibleActions,...
    possibleActions,possibleActions,possibleActions);

plotResults = true;                 % To plot or not the results at the end of the simulation
printInfo = false;                  % To print info after Bandits implementation (1) or not (0)
drawMap = false;                    % Variable for drawing the map
randomInitialConfiguration = true;  % Variable for assigning random channel/tx_power/cca at the beginning

% Dimensions of the 3D map
MaxX=10;
MaxY=5; 
MaxZ=10;
% Maximum range for a STA
MaxRangeX = 1;
MaxRangeY = 1;
MaxRangeZ = 1;
MaxRange = sqrt(3);

% Update modes of the exploration coefficient (epsilon-greedy)
UPDATE_MODE_FAST = 0;   % epsilon = initial_epsilon / t 
UPDATE_MODE_SLOW = 1;   % epsilon = epsilon / sqrt(t)

% Selected update mode (epsilon-greedy)
updateMode = UPDATE_MODE_SLOW;

% Path-loss model variables
PLd1=5;             % Path-loss factor
shadowing = 9.5;    % Shadowing factor
obstacles = 30;     % Obstacles factor
alfa = 4.4;         % Propagation model

COCHANNEL_INTERFERENCE = true;

%AoI
BW = 20e6;
num_info_bits = 100; %information bits length
frame_size = 200; %packet length
T = 240/BW;  %time slot

% Save constants into current folder
save('constants.mat');  