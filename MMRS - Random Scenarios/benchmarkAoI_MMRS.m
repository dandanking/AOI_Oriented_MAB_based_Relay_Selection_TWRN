function [tx_power_SR, tx_power_RD, instantaoiAfterAction_A, instantaoiAfterAction_B, time1, time2] = benchmarkAoI_MMRS(wlans,varargin)

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
    
    % Find the index of the initial action taken by each WLAN
    initialActionIxPerWlan = zeros(1, 2);
    for i = 1 : 2
        [~,indexCha] = find(channelActions == 1);
        [~,indexCca] = find(ccaActions == wlansAux(i).CCA);
        [~,indexTpc] = find(txPowerActions == wlansAux(i).TxPower);
        initialActionIxPerWlan(i) = indexes2val(indexCha, ...
            indexCca, indexTpc, size(channelActions,2), size(ccaActions,2));
    end

    SNRPerWlanA = zeros(2, K); 
    SNRPerWlanB = zeros(2, K); 
    actionIndexPerWlan = initialActionIxPerWlan;% Initialize the indexes of the taken action
    selectedArm = actionIndexPerWlan;
    timesArmHasBeenPlayed = zeros(2, K);  
    % channel_coefficient = 1/sqrt(2)*(randn(2,4) + j*randn(2,4));
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
    %The time of successful decoding
    time1 = zeros(1, totalIterations);
    time2 = zeros(1, totalIterations);
    % Record transmission power 
    tx_power_SR = zeros(2, totalIterations);
    tx_power_RD = zeros(2, totalIterations);
    %% ITERATE UNTIL CONVERGENCE OR MAXIMUM CONVERGENCE TIME       
    iteration = 1;
    rw = zeros(2,1);
    changeTime = totalIterations/10;
    while(iteration < totalIterations + 1) 
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
        for i =1:2
            if i == 1
                selectedArm(i) = select_action_MMRS(SNRPerWlanA, changeTime, iteration);
            else
                selectedArm(i) = select_action_MMRS(SNRPerWlanB, changeTime, iteration);
            end
            % Find channel and tx power of the current action
            [a, ~, c] = val2indexes(selectedArm(i), ...
                size(channelActions,2), size(ccaActions,2), size(txPowerActions,2)); 
            % Update WN configuration
            wlansAux(i).Channel = a;   
            wlansAux(i).TxPower = txPowerActions(c);
            timesArmHasBeenPlayed(i,selectedArm(i)) = timesArmHasBeenPlayed(i,selectedArm(i)) + 1;
            if i == 1
                SNRPerWlanA(1,selectedArm(i)) = pow2db(abs(channel_coefficient(1,a))^2 * db2pow(txPowerActions(c))) - NOISE_DBM;
                SNRPerWlanA(2,selectedArm(i)) = pow2db(abs(channel_coefficient(2,a))^2 * db2pow(1)) - NOISE_DBM;
            else
                SNRPerWlanB(1,selectedArm(i)) = pow2db(abs(channel_coefficient(2,a))^2 * db2pow(txPowerActions(c))) - NOISE_DBM;
                SNRPerWlanB(2,selectedArm(i)) = pow2db(abs(channel_coefficient(1,a))^2 * db2pow(1)) - NOISE_DBM;
            end
        end
        % Compute the instant aoi noticed after applying the action
        [CER_Wlan_A,CER_Wlan_B ] = compute_cer_from_sinr( wlansAux, channel_coefficient, SNRPerWlanA, SNRPerWlanB, selectedArm); %%%%%%%
        
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
         
        for i=1:2
            if wlansAux(1).Channel == wlansAux(2).Channel
                tx_power_SR(i,iteration) = wlansAux(i).TxPower;
                tx_power_RD(i,iteration) = 1;
            else
                tx_power_SR(i,iteration) = wlansAux(i).TxPower;
                tx_power_RD(i,iteration) = 2;
            end
        end
 
        
        % Increase the number of iterations
        iteration = iteration + 1;     
    end   
end

