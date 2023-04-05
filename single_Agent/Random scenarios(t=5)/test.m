clc
constants
% num_info_bits = 100; %information bits length
% frame_size = 200; %packet length
% for SNR=-5:2:10
%     sinr_dbm_SR = 10^(SNR/10);
%     fprintf('%d\n',sinr_dbm_SR);
%     %CER_SR_PerWlan1 = qfunc(log(2)*(log2(1+SNR)-num_info_bits/frame_size)/sqrt((1-1/(1+SNR)^2)/frame_size)); 
%     %CER_SR_PerWlan1 = qfunc(log(2)*(log2(1+sinr_dbm_SR)-num_info_bits/frame_size)/sqrt((1-1/(1+sinr_dbm_SR)^2)/frame_size)); 
%     CER_SR_PerWlan = qfunc((1/2*log2(1+10^(SNR/10))-num_info_bits/frame_size+log2(frame_size)/frame_size/2)/(log2(exp(1))*sqrt((1-1/(1+10^(SNR/10))^2)/frame_size/2))); %popovski complex AWGN
%     fprintf('%d %d\n',CER_SR_PerWlan1,CER_SR_PerWlan);
% end
% txPowerActions = [-2 0 2 4];
d1 = [20,30,40,50];
d2 = [50,40,30,20];
noise = -100;

for j=-15:5:5
    wlans(1).TxPower=j;
    fprintf('TxPower = %d\n', j);
    for i = 1:4
        d = d1(i);
        PL = PLd1 + 10 * alfa * log10(d) + shadowing / 2 + (d/10) .* obstacles / 2;
        powerStaFromAp_SR = wlans(1).TxPower - PL;  
        %fprintf('%.4f\n', powerStaFromAp_SR);
        noise = -100;
        % Compute the SINR per WLAN S->R                  
        sinr_dbm_SR = powerStaFromAp_SR - noise; %interference_plus_noise(i); % dBm
        fprintf('SNR = %.4f\n', sinr_dbm_SR);
    end
end
% 
% for m=-1:2
%     wlans(1).TxPower=m;
%     for n=-1:2
%         wlans(2).TxPower=n;
%         fprintf('S1.TxPower = %d S2.TxPower = %d\n', m,n);
%         for A_chn = 1:4
%             for B_chn = 1:4
%                 if A_chn ~=B_chn
%                     % Compute the received power on Relay      
%                     Ad = d1(A_chn);
%                     PL = PLd1 + 10 * alfa * log10(Ad) + shadowing / 2 + (Ad/10) .* obstacles / 2;
%                     powerStaFromAp_ASR = wlans(1).TxPower - PL;  
%                     Bd = d2(B_chn);
%                     PL1 = PLd1 + 10 * alfa * log10(Bd) + shadowing / 2 + (Bd/10) .* obstacles / 2;
%                     powerStaFromAp_BSR = wlans(2).TxPower - PL1;  
%                     % Compute the interference experienced per WLAN
%                     Ad = d2(A_chn);
%                     PL = PLd1 + 10 * alfa * log10(Ad) + shadowing / 2 + (Ad/10) .* obstacles / 2;
%                     interference_ASR = wlans(2).TxPower - PL;  
%                     Bd = d1(B_chn);
%                     PL1 = PLd1 + 10 * alfa * log10(Bd) + shadowing / 2 + (Bd/10) .* obstacles / 2;
%                     interference_BSR = wlans(1).TxPower - PL1; 
%                     interference_plus_noiseA = pow2db(db2pow(interference_ASR) + db2pow(noise));
%                     interference_plus_noiseB = pow2db(db2pow(interference_BSR) + db2pow(noise));
%                     % Compute the SINR per WLAN                   
%                     sinr_dbm_ASR = powerStaFromAp_ASR - interference_plus_noiseA; % dBm
%                     sinr_dbm_BSR = powerStaFromAp_BSR - interference_plus_noiseB; % dBm
%                     fprintf('S1—>R%d,SINR = %.4f,S2—>R%d, SINR = %.4f\n',A_chn, sinr_dbm_ASR,B_chn,sinr_dbm_BSR);
%                 end
%             end
%         end
%     end
% end

        