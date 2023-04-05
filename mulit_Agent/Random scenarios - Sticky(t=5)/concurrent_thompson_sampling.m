function [ aoiExperienced, timesArmHasBeenPlayed, tx_power_SR, tx_power_RD, instantaoiAfterAction_A, instantaoiAfterAction_B] = ...
    concurrent_thompson_sampling( wlans, varargin )
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
    
    actionIndexPerWlan = initialActionIxPerWlan;% Initialize the indexes of the taken action
    selectedArm = actionIndexPerWlan; % Initialize arm selection for each WLAN by using the initial action
    % Keep track of current and previous actions for getting the transitions probabilities
    currentAction = zeros(1, 2);
    previousAction = selectedArm;
    transitionsCounter = zeros(2, K^2);% Store the times a transition between actions is done in each WN
    timesArmHasBeenPlayed = zeros(2, K);% Store the times an action has been played in each WN
    cumulative_reward_per_action = zeros(2, K);%Cumulative reward per 10 action
    meanRewardPerWlan = zeros(2, K); % Initialize the mean reward obtained by each WLAN for each arm
    estimatedRewardPerWlan = ones(2, K)*(-10); % Initialize the mean reward obtained by each WLAN for each arm
    regretAfterAction = zeros(1, 2);% Initialize the regret experienced by each WLAN
    % random scenarios h
%             channel_coefficient = 1/sqrt(2)*(randn(2,4) + j*randn(2,4));
    channel_coefficient_1 = [-0.7700 + 0.3740i,2.0898 + 0.1620i,1.0928 + 0.8837i,-0.6938 + 0.7496i];
    channel_coefficient_1 = [channel_coefficient_1;channel_coefficient_1(:,end:-1:1)];
    channel_coefficient_2 = [-0.3474 + 0.7672i,0.1601 + 0.4126i,1.6944 + 0.5331i,0.1520 - 0.0351i];
    channel_coefficient_2 = [channel_coefficient_2;channel_coefficient_2(:,end:-1:1)];
    channel_coefficient_3 = [1.0823 + 0.4881i,1.2542 - 0.5718i,0.65171 - 0.0059i,0.2540 + 0.4475i];
    channel_coefficient_3 = [channel_coefficient_3;channel_coefficient_3(:,end:-1:1)];
    channel_coefficient_4 =[0.567552654409071 + 0.00405680991120112i,0.154971928876097 + 0.137202672687024i,-0.815021763384021 + 0.451885766888880i,-0.0434945684965565 - 0.536108102461996i];
    channel_coefficient_4 = [channel_coefficient_4;channel_coefficient_4(:,end:-1:1)];
    channel_coefficient_5 =[1.30183970901020 - 1.94773144431268i,-0.701870386900600 + 0.539789632888976i,1.06955581955932 - 0.352901658805943i,1.22223248675070 - 0.0790008655709236i];
    channel_coefficient_5 = [channel_coefficient_5;channel_coefficient_5(:,end:-1:1)];
    channel_coefficient_6 =[0.0448631617445665 + 1.15491166355308i,-0.554386269736583 + 0.452368105558358i,-1.17568088022720 - 0.611799528294491i,0.486767554025637 + 0.541670130016818i];
    channel_coefficient_6 = [channel_coefficient_6;channel_coefficient_6(:,end:-1:1)];
    channel_coefficient_7 =[-0.542710908931408 + 0.505394687160001i,-0.0332169590239790 + 0.0712739513031638i,0.271332305371790 + 0.227572618667519i,1.51552368636494 - 0.631815932473041i];
    channel_coefficient_7 = [channel_coefficient_7;channel_coefficient_7(:,end:-1:1)];
    channel_coefficient_8 =[0.313611582014330 - 0.440276162397075i,-0.635231041982387 + 0.222019484572058i,-0.929511288689500 - 0.352608572599453i,-0.695904611402084 + 0.118927302643564i];
    channel_coefficient_8 = [channel_coefficient_8;channel_coefficient_8(:,end:-1:1)];
    channel_coefficient_9 =[0.576411985181351 + 0.622515523681541i,-0.222725794209555 + 0.709618914422567i,1.76662936699024 + 0.878088263206310i,0.646950847446712 - 0.408020939342250i];
    channel_coefficient_9 = [channel_coefficient_9;channel_coefficient_9(:,end:-1:1)];
    channel_coefficient_10 =[1.59059931632655 - 0.167846199526370i,1.08212469520328 + 0.985924728168317i,0.0287675872643811 - 0.625717236267392i,0.169857173612801 + 1.12582436328979i];
    channel_coefficient_10 = [channel_coefficient_10;channel_coefficient_10(:,end:-1:1)];

    % Initialize the instant aoi per action
    instantaoiAfterAction_A = ones(1, totalIterations)*2;
    instantaoiAfterAction_B = ones(1, totalIterations)*2;
    %The time of successful decoding
    time1=zeros(1, totalIterations);
    time2=zeros(1, totalIterations);
    % Record transmission power 
    tx_power_SR=zeros(2, totalIterations);
    tx_power_RD=zeros(2, totalIterations);

    %% ITERATE UNTIL CONVERGENCE OR MAXIMUM CONVERGENCE TIME       
    iteration = 1;
    rw =ones(1,2).*(-10);
    while(iteration < totalIterations + 1) 
        % Assign turns to WLANs randomly 
        order = randperm(2);  
        if iteration <= 1000   
            channel_coefficient = channel_coefficient_1;
        elseif iteration <= 2000 
            channel_coefficient = channel_coefficient_2;
        elseif iteration <= 3000
            channel_coefficient = channel_coefficient_3;
        elseif iteration <= 4000
            channel_coefficient = channel_coefficient_4;
        elseif iteration <= 5000
            channel_coefficient = channel_coefficient_5;
        elseif iteration <= 6000
            channel_coefficient = channel_coefficient_6;
        elseif iteration <= 7000
            channel_coefficient = channel_coefficient_7;
        elseif iteration <= 8000
            channel_coefficient = channel_coefficient_8;
        elseif iteration <= 9000
            channel_coefficient = channel_coefficient_9;
        else
            channel_coefficient = channel_coefficient_10;
        end
        for i = 1:2
           if mod(iteration, 10) == 1 || iteration == 1
                %%%%%%%%%%%% Select an action according to the policy
                selectedArm(order(i)) = select_action_ts(estimatedRewardPerWlan(order(i), :), timesArmHasBeenPlayed(order(i), :));
                % Update the current action
                currentAction(order(i)) = selectedArm(order(i));
                % Find the index of the current and the previous action in allCombs
                ix = find(allCombs(:,1) == previousAction(order(i))...
                    & allCombs(:,2) == currentAction(order(i)));
                % Update the previous action
                previousAction(order(i)) = currentAction(order(i));  
                % Update the transitions counter
                transitionsCounter(order(i),ix) = transitionsCounter(order(i),ix) + 1;  
                % Find channel and tx power of the current action
                [a, ~, c] = val2indexes(selectedArm(order(i)), ...
                    size(channelActions,2), size(ccaActions,2), size(txPowerActions,2)); 
                % Update WN configuration
                wlansAux(order(i)).Channel = a;   
                wlansAux(order(i)).TxPower = txPowerActions(c);
           end
        end
        % Compute the instant aoi noticed after applying the action
        [CER_Wlan_A,CER_Wlan_B ] = compute_cer_from_sinr( wlansAux, NOISE_DBM, channel_coefficient); %%%%%%%
        
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
        end
        
        for wlan_i = 1:2
            % Update the times WN has selected the current action
            estimatedRewardPerWlan(wlan_i, selectedArm(wlan_i)) = ...
                (estimatedRewardPerWlan(wlan_i, selectedArm(wlan_i)) * timesArmHasBeenPlayed(wlan_i, selectedArm(wlan_i)) + rw(wlan_i)) /...
                (timesArmHasBeenPlayed(wlan_i, selectedArm(wlan_i)) + 2);
            cumulative_reward_per_action(wlan_i, selectedArm(wlan_i)) = cumulative_reward_per_action(wlan_i,selectedArm(wlan_i)) + rw(wlan_i);
            meanRewardPerWlan(wlan_i, selectedArm(wlan_i)) = ...
                cumulative_reward_per_action(wlan_i, selectedArm(wlan_i)) / timesArmHasBeenPlayed(wlan_i, selectedArm(wlan_i));
            timesArmHasBeenPlayed(wlan_i,selectedArm(wlan_i)) = ...
                    timesArmHasBeenPlayed(wlan_i,selectedArm(wlan_i)) + 1;
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