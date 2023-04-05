function [ aoiExperienced, timesArmHasBeenPlayed, tx_power_SR, tx_power_RD, instantaoiAfterAction_A, instantaoiAfterAction_B] = ...
    concurrent_egreedy_Sticky( wlans, initialEpsilon, varargin )
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
    load('channel_coefficient_(20_4)_1.mat');
    
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
    selectedAct = zeros(1, 2);
    for i = 1 : 2
        [~,indexCha] = find(channelActions == wlansAux(i).Channel);
        [~,indexCca] = find(ccaActions == wlansAux(i).CCA);
        [~,indexTpc] = find(txPowerActions == wlansAux(i).TxPower);
        initialActionIxPerWlan(i) = indexes2val(indexCha, ...
            indexCca, indexTpc, size(channelActions,2), size(ccaActions,2));
    end
    
    % Initialize the indexes of the taken action
    actionIndexPerWlan = initialActionIxPerWlan;                           
    % Initialize arm selection for each WLAN by using the initial action
    selectedArm = actionIndexPerWlan;     
    %selectedArm = [14 15]; 
    % Keep track of current and previous actions for getting the transitions probabilities
    currentAction = zeros(1, 2);
    previousAction = selectedArm;
    % Store the times a transition between actions is done in each WN
    transitionsCounter = zeros(2, K^2);    
    % Store the times an action has been played in each WN
    timesArmHasBeenPlayed = zeros(2, K^2);     
    % Initialize the mean reward obtained by each WLAN for each arm
    rewardPerArm_overall = zeros(2,K*K)*(-100);
    rewardPerArm_overall_AB = zeros(2,K*K)*(-100); 
    rewardPerArm = zeros(2, K).*(-100);   
    
    % random scenarios h
%             channel_coefficient = 1/sqrt(2)*(randn(2,4) + j*randn(2,4));
    channel_coefficient_1 = random_channel_coefficient(1:2,:);
    channel_coefficient_2 = random_channel_coefficient(3:4,:);
    channel_coefficient_3 = random_channel_coefficient(5:6,:);
    channel_coefficient_4 = random_channel_coefficient(7:8,:);
    channel_coefficient_5 = random_channel_coefficient(9:10,:);
    channel_coefficient_6 = random_channel_coefficient(11:12,:);
    channel_coefficient_7 = random_channel_coefficient(13:14,:);
    channel_coefficient_8 = random_channel_coefficient(15:16,:);
    channel_coefficient_9 = random_channel_coefficient(17:18,:);
    channel_coefficient_10 = random_channel_coefficient(19:20,:);
    
    % Initialize the instant aoi per action
    instantaoiAfterAction_A = ones(1, totalIterations)*2;
    instantaoiAfterAction_B = ones(1, totalIterations)*2;
    
    %Cumulative reward per 10 action
    cumulative_reward_per_action = zeros(2, K^2);

    %The time of successful decoding
    time1=zeros(1, totalIterations);
    time2=zeros(1, totalIterations);
    % Initialize epsilon
    epsilon = initialEpsilon; 
    % Record transmission power 
    tx_power_SR=zeros(2, totalIterations);
    tx_power_RD=zeros(2, totalIterations);

    %% ITERATE UNTIL CONVERGENCE OR MAXIMUM CONVERGENCE TIME       
    iteration = 1;
    rw =zeros(2,1);
    changeTime = totalIterations/10;
    while(iteration < totalIterations + 1) 
        % Assign turns to WLANs randomly 
        %order = randperm(nWlans);  
        if iteration <= changeTime   
%             channel_coefficient = 1/sqrt(2)*(randn(2,4) + j*randn(2,4));
            channel_coefficient = channel_coefficient_1;
        elseif iteration <= 2*changeTime 
            channel_coefficient = channel_coefficient_2;
        elseif iteration <= 3*changeTime
            channel_coefficient = channel_coefficient_3;
        elseif iteration <= 4*changeTime
            channel_coefficient = channel_coefficient_4;
        elseif iteration <= 5*changeTime
            channel_coefficient = channel_coefficient_5;
        elseif iteration <= 6*changeTime
            channel_coefficient = channel_coefficient_6;
        elseif iteration <= 7*changeTime
            channel_coefficient = channel_coefficient_7;
        elseif iteration <= 8*changeTime
            channel_coefficient = channel_coefficient_8;
        elseif iteration <= 9*changeTime
            channel_coefficient = channel_coefficient_9;
        else
            channel_coefficient = channel_coefficient_10;
        end
        for i = 1:2
           if wlansAux(i).Sticky <= 0
                if mod(iteration, 5) == 1 || iteration == 1
                    %%%%%%%%%%%% Select an action according to the policy
                    selectedArm(i) = select_action_egreedy(rewardPerArm_overall(i,:), epsilon); %%%%%改为很小
                end
           end
                % Update the current action
                currentAction(i) = selectedArm(i);
                % Find the index of the current and the previous action in allCombs
                ix = find(allCombs(:,1) == previousAction(i)...
                    & allCombs(:,2) == currentAction(i));
                % Update the previous action
                previousAction(i) = currentAction(i);  
                % Update the transitions counter
                transitionsCounter(i,ix) = transitionsCounter(i,ix) + 1;  
                % Find channel and tx power of the current action
                [selectedAct(1), ~,selectedAct(2)]=val2indexes(selectedArm(i),16,1,16);
                [a, ~, c] = val2indexes(selectedAct(i), ...
                    size(channelActions,2), size(ccaActions,2), size(txPowerActions,2)); 
                % Update WN configuration
                wlansAux(i).Channel = a;   
                wlansAux(i).TxPower = txPowerActions(c);
%            end
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
        if mod(iteration, 5) == 0 && iteration > 5 % from 20th
            for i=(iteration - 4):iteration
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
             rw(1) = - (sum1/10); 
            if rw(1) < benchmark_AOI
                wlansAux(1).Sticky = wlansAux(1).Sticky - 1; 
            else
                wlansAux(1).Sticky = 2;
            end
        end

        % Update the reward of WNB
        sum2 = 0;
        if mod(iteration, 5) == 0 && iteration > 5 % from 20th
            for i=(iteration - 4):iteration
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
             rw(2) = - (sum2/10); 
            if rw(2) < benchmark_AOI
                wlansAux(2).Sticky = wlansAux(2).Sticky - 1; 
            else
                wlansAux(2).Sticky = 2;
            end
        end
        
        for i = 1:2

            % Update the times WN has selected the current action
            timesArmHasBeenPlayed(i,selectedArm(i)) = ...
                    timesArmHasBeenPlayed(i,selectedArm(i)) + 1;
            rewardPerArm(i,selectedAct(i)) = rw(i);
            rewardPerArm_overall_AB(i, selectedArm(i))=rw(1)+rw(2);
            if ( CER_Wlan_B == 1) && (i == 1) % A收到B的历史轨迹
                rewardPerArm_overall(i,:) = rewardPerArm_overall_AB(i,:);
            end
            if ( CER_Wlan_A == 1) && (i == 2) % B收到A的历史轨迹
                rewardPerArm_overall(i,:) = rewardPerArm_overall_AB(i,:);
            end
            cumulative_reward_per_action(i, selectedArm(i)) = cumulative_reward_per_action(i,selectedArm(i)) + rw(i);
        end
        % Record transmission power 
        for i=1:2
            if wlansAux(1).Channel == wlansAux(2).Channel
                tx_power_SR(i,iteration) = wlansAux(i).TxPower;
                tx_power_RD(i,iteration) = 1;
            else
                tx_power_SR(i,iteration) = wlansAux(i).TxPower;
                tx_power_RD(i,iteration) = 2;
            end
        end

        % Store the throughput and the regret at the end of the iteration for statistics
        aoiExperienced(iteration,:) = -rw;

        % Update the exploration coefficient according to the inputted mode
        if updateMode == UPDATE_MODE_FAST
            epsilon = initialEpsilon / iteration;    
        elseif updateMode == UPDATE_MODE_SLOW
            epsilon = initialEpsilon / sqrt(iteration);   
        elseif updateMode == UPDATE_MODE_CONSTANT
            epsilon = 0.02;
        else
            disp(['updateModeEpsilon = ' num2str(updateModeEpsilon) ' does not exist!'])
        end        
        % Increase the number of iterations
        iteration = iteration + 1;     
    end   
    for i = 1:2
        for action_ix = 1 : K
            meanRewardPerAction(i,action_ix) =  ...
                cumulative_reward_per_action(i,action_ix)...
                / timesArmHasBeenPlayed(i,action_ix);
        end    
    end

end