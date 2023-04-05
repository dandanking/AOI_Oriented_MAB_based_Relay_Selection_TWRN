function [ aoiExperienced, timesArmHasBeenPlayed, tx_power_SR, tx_power_RD, instantaoiAfterAction_A, instantaoiAfterAction_B] = ...
    concurrent_exp3( wlans, gamma, initialEta, varargin )
% EGREEDY - Given a WN, applies e-greedy to maximize the experienced throughput
%
%   OUTPUT: 
%       * aoiExperiencedByWlan - throughput experienced by each WLAN
%         for each of the iterations done
%       * timesArmHasBeenPlayed - times each action has been played
%   INPUT: 
%       * wlan - wlan object containing information about all the WLANs
%       * initialEpsilon - initial exploration coefficient

    load('constants.mat')
    
    try
        if size(varargin, 2) == 3
            % Update possible actions
            nChannels = varargin{1};
            channelActions = 1 : nChannels;
            ccaActions = varargin{2};
            txPowerActions = varargin{3};
            % Each state represents an [i,j,k] combination for indexes on "channels", "cca" and "tx_power"
            possibleActions = 1:(size(channelActions, 2) * ...
                size(ccaActions, 2) * size(txPowerActions, 2));
            K = size(possibleActions,2);   % Total number of actions
            allCombs = allcomb(1:K, 1:K);    
        end
    catch
        if size(varargin, 2) ~= 3
            disp('Wrong number of input arguments')
        end
    end
    
    %% INITIALIZE ALGORITHM
    % Use a copy of wlan to make operations
    wlansAux = wlans;
    nWlans = size(wlansAux, 2);
    
    % Find the index of the initial action taken by each WLAN
    initialActionIxPerWlan = zeros(1, 2);
    for i = 1 : 2
        [~,indexCha] = find(channelActions == wlansAux(i).Channel);
        [~,indexCca] = find(ccaActions == wlansAux(i).CCA);
        [~,indexTpc] = find(txPowerActions == wlansAux(i).TxPower);
        initialActionIxPerWlan(i) = indexes2val(indexCha, ...
            indexCca, indexTpc, size(channelActions,2), size(ccaActions,2));
    end
    
    actionIndexPerWlan = initialActionIxPerWlan; % Initialize the indexes of the taken action
    selectedArm = actionIndexPerWlan;   % Initialize arm selection for each WLAN by using the initial action
    
    weightsPerArm = ones(2, K);    % Initialize weight to 1 for each action
    previousAction = selectedArm;       % Keep track of current and previous actions for getting the transitions probabilities
    timesArmHasBeenPlayed = zeros(2, K);% Store the times an action has been played in each WN
    transitionsCounter = zeros(2, K^2); % Store the times a transition between actions is done in each WN   
    
    armsProbabilities = (1/K)*ones(2, K); % Initialize arms probabilities
    estimated_reward = zeros(1, 2);       % Initialize the estimated reward for each WN

    % Initialize the instant aoi per action
    instantaoiAfterAction_A = ones(1, totalIterations)*2;
    instantaoiAfterAction_B = ones(1, totalIterations)*2;
    %The time of successful decoding
    time1=zeros(1, totalIterations);
    time2=zeros(1, totalIterations);
    % Record transmission power 
    tx_power_SR=zeros(2, totalIterations);
    tx_power_RD=zeros(2, totalIterations);

    % Initialize the learning rate
    eta = initialEta;
    previousEta = eta;    
    %% ITERATE UNTIL CONVERGENCE OR MAXIMUM CONVERGENCE TIME       
    iteration = 1;
    rw =zeros(1,2);
    while(iteration < totalIterations + 1) 
        % Assign turns to WLANs randomly 
        %order = randperm(nWlans);  
        for i = 1:2
           if mod(iteration, 10) == 1 || iteration == 1
                % Update arms probabilities according to weights      
                for k = 1 : K
                    armsProbabilities(i, k) = (1 - gamma) * (weightsPerArm(i, k) / sum(weightsPerArm(i, :))) + (gamma / K);                
                    % To avoid errors in execution time
                    if isnan(armsProbabilities(i, k))
                        armsProbabilities(i, k) = 1/K;
                    end
                end 
                % Draw an action according to probabilites distribution
                selectedArm(i) = randsample(1:K, 1, true, armsProbabilities(i,:));  
                % Find the index of the current and the previous action in allCombs
                ix = find(allCombs(:,1) == previousAction(i) & allCombs(:,2) == selectedArm(i));
                % Update the previous action
                previousAction(i) = selectedArm(i);  
                % Update the transitions counter
                transitionsCounter(i,ix) = transitionsCounter(i,ix) + 1;  
                % Update the times WN has selected the current action
                timesArmHasBeenPlayed(i,selectedArm(i)) = timesArmHasBeenPlayed(i,selectedArm(i)) + 1;                
                % Find channel and tx power of the current action
                [a, ~, c] = val2indexes(selectedArm(i), size(channelActions,2), size(ccaActions,2), size(txPowerActions,2)); 
                % Update WN configuration
                wlansAux(i).Channel = a;   
                wlansAux(i).TxPower = txPowerActions(c);
           end
        end
        % Compute the instant aoi noticed after applying the action
        [CER_Wlan_A,CER_Wlan_B ] = compute_cer_from_sinr( wlansAux, NOISE_DBM, iteration); %%%%%%%
        
        %UserA
        if ( CER_Wlan_A == 0) && (iteration >= 2) %send failure
            instantaoiAfterAction_A(iteration) = instantaoiAfterAction_A(iteration-1) + 2;
        else
            instantaoiAfterAction_A(iteration) = instantaoiAfterAction_A(iteration) + 2;
            time1(iteration) = 1;
        end
        %UserB
        if ( CER_Wlan_B == 0) && (iteration >= 2) %send failure
            instantaoiAfterAction_B(iteration) = instantaoiAfterAction_B(iteration-1) + 2;
        else
            instantaoiAfterAction_B(iteration) = instantaoiAfterAction_B(iteration) + 2;
            time2(iteration) = 1;
        end
        
        % Update the reward of WNA
        sum1 = 0;
        if mod(iteration, 10) == 0 && iteration > 10 % from 20th
            for i=(iteration - 9):iteration
                if time1(i)==1 && time1(i-1)==1
                    sum1 = sum1 + (instantaoiAfterAction_A(i)+2)*2/2;
                end
                if time1(i) ==1 && time1(i-1) ~=1
                    sum1 = sum1 + (instantaoiAfterAction_A(i)+instantaoiAfterAction_A(i-1))*2/2;
                end
                if time1(i) ~= 1 && time1(i-1)==1
                    sum1 = sum1 + (instantaoiAfterAction_A(i)+2)*2/2;
                end
                 if time1(i) ~= 1 && time1(i-1) ~= 1
                    sum1 = sum1 + (instantaoiAfterAction_A(i)+instantaoiAfterAction_A(i-1))*2/2;
                end
            end
             rw(1) = - (sum1/20); 
        end

        % Update the reward of WNB
        sum2 = 0;
        if mod(iteration, 10) == 0 && iteration > 10 % from 20th
            for i=(iteration - 9):iteration
                if time2(i)==1 && time2(i-1)==1
                    sum2 = sum2 + (instantaoiAfterAction_B(i)+2)*2/2;
                end
                if time2(i) ==1 && time2(i-1) ~=1
                    sum2 = sum2 + (instantaoiAfterAction_B(i)+instantaoiAfterAction_B(i-1))*2/2;
                end
                if time2(i) ~= 1 && time2(i-1)==1
                    sum2 = sum2 + (instantaoiAfterAction_B(i)+2)*2/2;
                end
                 if time2(i) ~= 1 && time2(i-1) ~= 1
                    sum2 = sum2 + (instantaoiAfterAction_B(i)+instantaoiAfterAction_B(i-1))*2/2;
                end
            end
             rw(2) = - (sum2/20); 
        
            for wlan_i = 1 : 2           
                % Update the estimated reward
                estimated_reward(wlan_i) = (rw(wlan_i) / armsProbabilities(wlan_i, selectedArm(wlan_i)));                    
                % Update the weights of eah action
                for k = 1 : K
                    if eta == 0 && previousEta == 0
                        weightsPerArm(wlan_i, k) = weightsPerArm(wlan_i, k)^0 * ...
                            exp((eta * estimated_reward(wlan_i)));           
                    else
                        if k == selectedArm(wlan_i) 
                            weightsPerArm(wlan_i, k) = weightsPerArm(wlan_i, k)^...
                                (eta / previousEta) * exp((eta * estimated_reward(wlan_i)));
                        else
                            weightsPerArm(wlan_i, k) = ...
                                weightsPerArm(wlan_i, k)^(eta / previousEta);
                        end
                     end
                     weightsPerArm(wlan_i, k) = max( weightsPerArm(wlan_i, k), 1e-6 );
                end
            end 
        end
        % Record transmission power 
        for i=1:2
            if wlansAux(1).Channel == wlansAux(2).Channel
                tx_power_SR(i,iteration) = wlansAux(i).TxPower;
                tx_power_RD(i,iteration) = 2;
            else
                tx_power_SR(i,iteration) = wlansAux(i).TxPower;
                tx_power_RD(i,iteration) = 4;
            end
        end

        % Store the throughput and the regret at the end of the iteration for statistics
        aoiExperienced(iteration, :) = -rw; 

        % Update the learning rate according to the "update mode"
        previousEta = eta;
        if updateMode == UPDATE_MODE_FAST
            eta = initialEta / iteration;    
        elseif updateMode == UPDATE_MODE_SLOW
            eta = initialEta / sqrt(iteration);   
        else
            % eta remains constant	
        end     
        % Increase the number of iterations
        iteration = iteration + 1;     
    end   
    %% PRINT INFORMATION REGARDING ACTION SELECTION
    if printInfo    
        % Print the preferred action per wlan
        for i= 1 : nWlans      
            timesArmHasBeenPlayed(i, :)/totalIterations
            a = transitionsCounter(i,:);
            % Max value
            [val1, ix1] = max(a);
            [ch1_1, ~, x] = val2indexes(possibleActions(allCombs(ix1,1)), size(channelActions,2), size(ccaActions,2), size(txPowerActions,2)); 
            tpc1_1 = txPowerActions(x);
            [ch1_2, ~, x] = val2indexes(possibleActions(allCombs(ix1,2)), size(channelActions,2), size(ccaActions,2), size(txPowerActions,2)); 
            tpc1_2 = txPowerActions(x);
            % Second max value
            [val2, ix2] = max(a(a<max(a)));
            [ch2_1, ~, x] = val2indexes(possibleActions(allCombs(ix2,1)), size(channelActions,2), size(ccaActions,2), size(txPowerActions,2)); 
            tpc2_1 = txPowerActions(x);
            [ch2_2, ~, x] = val2indexes(possibleActions(allCombs(ix2,2)), size(channelActions,2), size(ccaActions,2), size(txPowerActions,2)); 
            tpc2_2 = txPowerActions(x);
            % Third max value
            [val3, ix3] = max(a(a<max(a(a<max(a)))));
            [ch3_1, ~, x] = val2indexes(possibleActions(allCombs(ix3,1)), size(channelActions,2), size(ccaActions,2), size(txPowerActions,2)); 
            tpc3_1 = txPowerActions(x);
            [ch3_2, ~, x] = val2indexes(possibleActions(allCombs(ix3,2)), size(channelActions,2), size(ccaActions,2), size(txPowerActions,2)); 
            tpc3_2 = txPowerActions(x);   

            disp(['Probability of going from ' num2str(allCombs(ix1,1)) ' (ch=' num2str(ch1_1) '/tpc=' num2str(tpc1_1) ')' ...
                ' to ' num2str(allCombs(ix1,2)) ' (ch=' num2str(ch1_2) '/tpc=' num2str(tpc1_2) ')' ...
                ' = ' num2str(val1/totalIterations)])

            disp(['Probability of going from ' num2str(allCombs(ix2,1)) ' (ch=' num2str(ch2_1) '/tpc=' num2str(tpc2_1) ')' ...
                ' to ' num2str(allCombs(ix2,2)) ' (ch=' num2str(ch2_2) '/tpc=' num2str(tpc2_2) ')' ...
                ' = ' num2str(val2/totalIterations)])

            disp(['Probability of going from ' num2str(allCombs(ix3,1)) ' (ch=' num2str(ch3_1) '/tpc=' num2str(tpc3_1) ')' ...
                ' to ' num2str(allCombs(ix3,2)) ' (ch=' num2str(ch3_2) '/tpc=' num2str(tpc3_2) ')' ...
                ' = ' num2str(val3/totalIterations)])

        end        
    end
    
end