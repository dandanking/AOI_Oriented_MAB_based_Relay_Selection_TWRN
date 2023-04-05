function [snr_db_ASR,snr_db_BSR]  = compute_rpower_per_wlan(wlans,noise)
% Interferences - Returns the interferences power received at each WLAN
%   OUTPUT:
%       * intMat: 1xN array (N is the number of WLANs) with the
%       receive power noticed on each AP in mW
%   INPUT:
%       * wlans: contains information of each WLAN in the map. For instance,
%       wlans(1) corresponds to the first one, so that it has unique
%       parameters (x,y,z,BW,CCA,etc.).

    wlans_aux = wlans;
   
    load('constants.mat')
    
    
    powerStaFromAp_ASR = zeros(1, size(wlans_aux, 2)); 
    powerStaFromAp_BSR = zeros(1, size(wlans_aux, 2)); 
    snr_db_ASR = zeros(1, size(wlans_aux, 2));
    snr_db_BSR = zeros(1, size(wlans_aux, 2));
    d1 = [17,18,19,20];
    d2 = [20,19,18,17];
    % Compute the received power on Relay  
    for i=1:4
        Ad = d1(i);
        PL = PLd1 + 10 * alfa * log10(Ad) + shadowing / 2 + (Ad/10) .* obstacles / 2;
        powerStaFromAp_ASR(i) = wlans(1).TxPower - PL;  
        Bd = d2(i);
        PL1 = PLd1 + 10 * alfa * log10(Bd) + shadowing / 2 + (Bd/10) .* obstacles / 2;
        powerStaFromAp_BSR(i) = wlans(2).TxPower - PL1;  
        snr_db_ASR(i) = powerStaFromAp_ASR(i) - noise; 
        snr_db_BSR(i) = powerStaFromAp_BSR(i) - noise; 
    end
%     Ad = d1(B_chn);
%     PL = PLd1 + 10 * alfa * log10(Ad) + shadowing / 2 + (Ad/10) .* obstacles / 2;
%     interference_ASR = wlans(2).TxPower - PL;  
%     Bd = d2(A_chn);
%     PL1 = PLd1 + 10 * alfa * log10(Bd) + shadowing / 2 + (Bd/10) .* obstacles / 2;
%     interference_BSR = wlans(1).TxPower - PL1; 
% 
%     interference_plus_noiseA = pow2db(interference_ASR + db2pow(noise));
%     interference_plus_noiseB = pow2db(interference_BSR + db2pow(noise));

    % Compute the SNR per WLAN                   

%     for i = 1 : num_wlans
%         for j = 1 : num_wlans
%             if i ~= j 
%                 if COCHANNEL_INTERFERENCE
%                     adjacent_interference = 20 * (abs(wlans_aux(i).Channel - wlans_aux(j).Channel));
%                     if wlans_aux(i).Channel ~= wlans_aux(j).Channel
%                         interference_per_wlan_mw(i) = interference_per_wlan_mw(i) + db2pow(powerMatrix(i,j) - adjacent_interference);
%                     else
%                         interference_per_wlan_mw(i) = interference_per_wlan_mw(i) + db2pow(powerMatrix(i,j));
%                     end  
%                 else
%                     if wlans_aux(i).Channel == wlans_aux(j).Channel
%                         interference_per_wlan_mw(i) = interference_per_wlan_mw(i) + db2pow(powerMatrix(i,j));
%                     end
%                 end
%             end
%         end
%     end

end