function [CER_Wlan_A,CER_Wlan_B] = compute_cer_from_sinr( wlans, channel_coefficient, SNRPerWlanA, SNRPerWlanB, selectedArm)
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
    
    snr_db_AR = SNRPerWlanA(1,selectedArm(1)); 
    snr_db_RA = SNRPerWlanA(2,selectedArm(1));
    snr_db_BR = SNRPerWlanB(1,selectedArm(2));
    snr_db_RB = SNRPerWlanB(2,selectedArm(2));
    
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