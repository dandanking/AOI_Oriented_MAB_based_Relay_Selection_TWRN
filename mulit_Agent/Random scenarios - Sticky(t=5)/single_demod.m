function [ ABits, p_diff_A] ...
    = single_demod(rx, ha, modulation, pkt_length, snr )
logp_A = zeros(pkt_length,1);
p_diff_A = zeros(pkt_length,1);
pA_0 = zeros(pkt_length,1);
pA_1 = zeros(pkt_length,1);
%compute probability
for ii=1:length(rx)
    snr_lin=10^(-snr/10);
    if strcmp(modulation,'BPSK') == 1
        data_pair = [
            1;
            -1];   
    end    
    d = zeros(length(data_pair),1);
    for pair_index = 1 : length(data_pair)
        d(pair_index) = norm( ha(ii) * data_pair(pair_index, 1)  ...
            - rx(ii) )^2;
    end    
    if strcmp(modulation,'BPSK') == 1
        %packet A
        p0_A = exp(-d(2)^2/snr_lin);
        p1_A = exp(-d(1)^2/snr_lin);
        p_diff_A(ii) = p0_A-p1_A;
        logp_A(ii) = log(p0_A-p1_A);
        pA_0(ii) = p0_A;
        pA_1(ii) = p1_A;
    end    
    ABits = -1*reshape([pA_0,pA_1]',1,length(pA_0)*2);
end
end