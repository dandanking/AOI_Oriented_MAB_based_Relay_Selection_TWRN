function [CER_Wlan_A,CER_Wlan_B] = compute_cer_from_sinr( wlans, NOISE_DBM, channel_coefficient)
% Computes the cer in wlan 
% OUTPUT:
%   * CER_Wlan 
% INPUT:
%   * wlan - object containing all the WLANs information 

    constants
    % 变量初始化
    
    Pd_SR = zeros(2, 1);
    A_chn = wlans(1).Channel;
    B_chn = wlans(2).Channel;
    A_txPow = wlans(1).TxPower;
    B_txPow = wlans(2).TxPower;
    
    snr_db_AR = pow2db(abs(channel_coefficient(1,A_chn))^2 * db2pow(A_txPow)) - NOISE_DBM;
    snr_db_RA = pow2db(abs(channel_coefficient(2,A_chn))^2 * db2pow(1)) - NOISE_DBM;
    snr_db_BR = pow2db(abs(channel_coefficient(2,B_chn))^2 * db2pow(B_txPow)) - NOISE_DBM;
    snr_db_RB = pow2db(abs(channel_coefficient(1,B_chn))^2 * db2pow(1)) - NOISE_DBM;
    
    %% Select the same relay
    if A_chn == B_chn
        % 第一跳
        [okX] = PNC_phy_sim(snr_db_AR, snr_db_BR, channel_coefficient(1,A_chn), channel_coefficient(2,B_chn));    %%%%%%%%%%%%%
        
        % 第二跳
        if okX == 1
            [ okA] = single_phy_sim(snr_db_RA, channel_coefficient(2,A_chn)); %%%%%%%%%
            [ okB] = single_phy_sim(snr_db_RB, channel_coefficient(1,B_chn)); %%%%%%%%%
            CER_Wlan_A = okA;
            CER_Wlan_B = okB;
        else
            CER_Wlan_A = 0;
            CER_Wlan_B = 0;
        end
   %% Select different relay
    else
        % first-hop
        [ okA] = single_phy_sim(snr_db_AR, channel_coefficient(1,A_chn));
        Pd_SR(1) = okA;
        [ okA] = single_phy_sim(snr_db_BR, channel_coefficient(2,B_chn));
        Pd_SR(2) = okA;
        
        % second-hop
        if (Pd_SR(1) == 1)
            [okXA] = single_phy_sim(snr_db_RA, channel_coefficient(2,A_chn));
            CER_Wlan_A = okXA;
        else
            CER_Wlan_A = 0;
        end
        if (Pd_SR(2) == 1)
            [okXB] = single_phy_sim(snr_db_RB, channel_coefficient(1,B_chn));
            CER_Wlan_B = okXB;
        else
            CER_Wlan_B = 0;
        end
   end    
end