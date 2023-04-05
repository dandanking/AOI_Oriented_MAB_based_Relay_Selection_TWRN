function arm = select_action_MMRS(SNRPerWlan, changeTime, iter)
%   OUTPUT:
%       * arm - which represents the configuration composed by [channel,CCA,TPC]
%   INPUT:
%       * channel_coefficient_per_configuration - channel coefficient
%       between S and R(R and S)

    indexes=[];  % empty array to include configurations with the same value    
    
    if mod(iter,changeTime) > 16  || mod(iter,changeTime) == 0 % max-min

        min_SNR = min(SNRPerWlan,[],1);   %选择每列最小值
        [val,~] = max(min_SNR);           %选择最大值
        
        % Break ties randomly
        if sum(min_SNR==val)>1
            if val ~= Inf
                indexes = find(min_SNR==val);
                arm = randsample(indexes,1);
            else
                arm = randsample(1:size(min_SNR,2),1);
            end
            
        % Select arm with maximum min_SNR
        else
            [~,arm] = max(min_SNR);
        end
        
    else         % 遍历     

        arm = mod(iter,changeTime);
        
    end
end