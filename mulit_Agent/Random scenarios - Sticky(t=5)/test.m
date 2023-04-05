clc
% constants
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
% d1 = [20,30,40,50];
% d2 = [50,40,30,20];
% noise = -100;
% 
% for j=-15:5:5
%     wlans(1).TxPower=j;
%     fprintf('TxPower = %d\n', j);
%     for i = 1:4
%         d = d1(i);
%         PL = PLd1 + 10 * alfa * log10(d) + shadowing / 2 + (d/10) .* obstacles / 2;
%         powerStaFromAp_SR = wlans(1).TxPower - PL;  
%         %fprintf('%.4f\n', powerStaFromAp_SR);
%         noise = -100;
%         % Compute the SINR per WLAN S->R                  
%         sinr_dbm_SR = powerStaFromAp_SR - noise; %interference_plus_noise(i); % dBm
%         fprintf('SNR = %.4f\n', sinr_dbm_SR);
%     end
% end
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

    % random scenarios h
    channel_coefficient = 1/sqrt(2)*(randn(20,4) + j*randn(20,4));
    channel_coefficient_1 = [0.246770250942284 - 1.17628456790003i,0.238926054375777 - 0.0432085142899853i,1.30032162558941 - 0.540959325916110i,-0.695874443963344 + 0.586870337848432i;
       -0.403681003070308 - 0.107067656595967i,0.487657349300139 - 1.17511913521600i,0.647756149260672 - 1.47429936714327i,0.427188174791745 - 0.456940402028544i];
    channel_coefficient_2 = [0.263608023522948 - 0.0867375770981208i,0.884978019998268 + 1.62334247868670i,0.00883339177816810 + 0.467768834673257i,-0.725476834579955 - 0.0299151960875711i;
        -0.415627280440092 - 0.0288098286467097i,-1.15476963402491 + 0.397260252259251i,-0.783429220682272 - 0.867793374062022i,0.951452306124265 - 0.134416152125293i];
    channel_coefficient_3 = [0.328660594969264 - 1.82791813319194i,-0.953953422800092 - 0.567421217495743i,1.77120094623109 - 0.170943341746866i,-0.0596158023972971 - 0.195811615340043i;
       0.440792286484106 - 0.492955408269866i,-0.0102913615856595 + 0.450592331608364i,0.951864342707867 + 0.456988723264963i,0.830526237489649 - 0.438311538225874i];
    channel_coefficient_4 =[-0.232988530071136 - 0.662314474780249i,-0.474873838095093 - 0.100221207362854i,-0.961066408062129 - 0.321480106422286i,0.389569903877464 + 0.457911709978447i;
        0.332498406559051 + 0.783398027113733i,0.481573790432486 - 0.748800685768488i,-0.170698553099505 + 0.142305523981088i,-0.524133650369548 - 0.0312719173042769i];
    channel_coefficient_5 =[-1.07387307913638 + 0.991808700162391i,0.0634356579427836 - 0.234388589023115i,-0.126139117311536 - 1.75182219193877i,0.525871809037861 - 1.43123069066487i;
        -0.285854441615477 - 0.560271048270150i,0.476821627242443 + 0.592768147502002i,-1.08176949082250 - 0.401717272718555i,-0.828714561866430 - 0.506736114913743i];
    channel_coefficient_6 =[-0.472599109798949 + 0.236632602703155i,0.847019806190732 - 0.477751511520269i,-1.31724650419445 + 1.15706594954012i,0.717567467461353 + 0.155581843418036i;
        -0.0105148597493184 + 0.336063872535342i,0.689834777672800 + 0.591768389499996i,0.705692451975618 - 0.115782722323992i,-0.408163297743940 - 0.489308024128660i];
    channel_coefficient_7 =[-0.0964822888483367 + 1.11243223484887i,-0.0583880385841197 - 0.0788988151568244i,-2.04687019844307 + 1.30947951861100i,-0.279958685218910 - 0.562942867634598i;
        0.979038479731059 + 0.220651534148718i,-0.601377973900525 - 0.956290106784870i,-1.38134473920359 + 0.950277861089146i,-0.0766291432979004 + 0.199510194674669i];
    channel_coefficient_8 =[-0.0436591086750099 + 0.104956712925688i,0.713084490342131 + 1.43204062830614i,0.719276279225745 + 0.407375427508203i,-0.748540105856548 + 0.222871696009562i;
        0.535031255313341 - 0.713752676615435i,0.186938341908351 - 0.0915617066976199i,0.0379215657692034 + 0.177514191711558i,0.102106604583449 + 0.264390986838590i];
    channel_coefficient_9 =[-1.27251316579682 - 0.741943194025465i,0.820744953829873 - 0.384257744239362i,0.284498311124554 - 0.542091188706067i,-1.40557826748126 + 0.116242793403601i;
        -1.25087477974509 + 0.282036753445749i,0.989566706742766 - 0.613827325170720i,-0.444261723613850 - 0.0728905172626306i,0.740419572078394 + 0.303173419372990i];
    channel_coefficient_10 =[-1.39118512742771 + 0.0469235379147789i,0.397934122812516 - 0.395584537723068i,1.60101398184773 - 0.0204186687592850i,-0.0464359762470070 + 0.785236804442400i;
        0.385715089712530 - 0.0231152905099388i,-0.249609203903836 + 1.06066208007810i,-0.0314573420626216 - 0.0511996971150158i,0.599435687169190 - 0.0384724966940222i];
