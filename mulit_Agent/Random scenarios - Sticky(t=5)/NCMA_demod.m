function [ ABits, BBits, XBits,p_diff_A,p_diff_B,p_diff_xor ] ...
    = NCMA_demod(rx, ha, hb, modulation, pkt_length, snr )

logp_A = zeros(pkt_length,1);
p_diff_A = zeros(pkt_length,1);

logp_B = zeros(pkt_length,1);
p_diff_B = zeros(pkt_length,1);

logp_xor = zeros(pkt_length,1);
p_diff_xor = zeros(pkt_length,1);

pA_0 = zeros(pkt_length,1);
pA_1 = zeros(pkt_length,1);
pB_0 = zeros(pkt_length,1);
pB_1 = zeros(pkt_length,1);
pX_0 = zeros(pkt_length,1);
pX_1 = zeros(pkt_length,1);

%compute probability
for ii=1:length(rx)
    snr_lin=10^(-snr/10);
    if strcmp(modulation,'BPSK') == 1
        data_pair = [1 1;
            1 -1;
            -1 1;
            -1 -1];
    elseif strcmp(modulation,'QPSK') == 1
        data_pair= [
            1+j  1+j;
            1+j  -1+j;
            1+j  1-j;
            1+j  -1-j;
            
            -1+j  1+j;
            -1+j  -1+j;
            -1+j  1-j;
            -1+j  -1-j;
            
            1-j  1+j;
            1-j  -1+j;
            1-j  1-j;
            1-j  -1-j;
            
            -1-j  1+j;
            -1-j  -1+j;
            -1-j  1-j;
            -1-j  -1-j]/sqrt(2);
        
    end
    
    d = zeros(length(data_pair),1);
    for pair_index = 1 : length(data_pair)
        d(pair_index) = norm( ha(ii) * data_pair(pair_index, 1)  ...
            + hb(ii) * data_pair(pair_index,2)  ...
            - rx(ii) )^2;

    end
    
    if strcmp(modulation,'BPSK') == 1
        %packet A
        p0_A = exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin);
        p1_A = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin);
        p_diff_A(ii) = p0_A-p1_A;
        logp_A(ii) = log(p0_A-p1_A);
        pA_0(ii) = p0_A;
        pA_1(ii) = p1_A;
        
        
        %packet B
        p0_B = exp(-d(2)^2/snr_lin)+exp(-d(4)^2/snr_lin);
        p1_B = exp(-d(1)^2/snr_lin)+exp(-d(3)^2/snr_lin);
        p_diff_B(ii) = p0_B-p1_B;
        logp_B(ii) = log(p0_B-p1_B);
        pB_0(ii) = p0_B;
        pB_1(ii) = p1_B;
        
        %packet xor
        p0_xor = exp(-d(1)^2/snr_lin)+exp(-d(4)^2/snr_lin);
        p1_xor = exp(-d(2)^2/snr_lin)+exp(-d(3)^2/snr_lin);
        p_diff_xor(ii) = p0_B-p1_B;
        logp_xor(ii) = log(p0_xor-p1_xor);
        pX_0(ii) = p0_xor;
        pX_1(ii) = p1_xor;
        
        
    elseif strcmp(modulation,'QPSK') == 1
        
        %packet A
        p0_A_I = exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin)+...
            exp(-d(13)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(15)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_A_I = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin)+...
            exp(-d(9)^2/snr_lin)+exp(-d(10)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(12)^2/snr_lin);
        p_diff_A(2*ii-1) = p0_A_I-p1_A_I;
        logp_A(2*ii-1) = log(p0_A_I-p1_A_I);
        pA_0(2*ii-1) = p0_A_I;
        pA_1(2*ii-1) = p1_A_I;
        
        p0_A_Q = exp(-d(9)^2/snr_lin)+exp(-d(10)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(12)^2/snr_lin)+...
            exp(-d(13)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(15)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_A_Q = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin)+...
            exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin);
        p_diff_A(2*ii) = p0_A_Q-p1_A_Q;
        logp_A(2*ii) = log(p0_A_Q-p1_A_Q);
        pA_0(2*ii) = p0_A_Q;
        pA_1(2*ii) = p1_A_Q;
        
        
        %packet B
        p0_B_I = exp(-d(2)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(8)^2/snr_lin)+...
            exp(-d(10)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_B_I = exp(-d(1)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(7)^2/snr_lin)+...
            exp(-d(9)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(15)^2/snr_lin);
        p_diff_B(2*ii-1) = p0_B_I-p1_B_I;
        logp_B(2*ii-1) = log(p0_B_I-p1_B_I);
        pB_0(2*ii-1) = p0_B_I;
        pB_1(2*ii-1) = p1_B_I;
        
        p0_B_Q = exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin)+...
            exp(-d(11)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(15)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_B_Q = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin)+...
            exp(-d(9)^2/snr_lin)+exp(-d(10)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(14)^2/snr_lin);
        p_diff_B(2*ii) = p0_B_Q-p1_B_Q;
        logp_B(2*ii) = log(p0_B_Q-p1_B_Q);
        pB_0(2*ii) = p0_B_Q;
        pB_1(2*ii) = p1_B_Q;
        
        %packet xor
        p0_xor_I = exp(-d(1)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(8)^2/snr_lin)+...
            exp(-d(9)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xor_I = exp(-d(2)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(7)^2/snr_lin)+...
            exp(-d(10)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(15)^2/snr_lin);
        p_diff_xor(2*ii-1) = p0_xor_I-p1_xor_I;
        logp_xor(2*ii-1) = log(p0_xor_I-p1_xor_I);
        pX_0(2*ii-1) = p0_xor_I;
        pX_1(2*ii-1) = p1_xor_I;
        
        p0_xor_Q = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin)+...
            exp(-d(11)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(15)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xor_Q = exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin)+...
            exp(-d(9)^2/snr_lin)+exp(-d(10)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(14)^2/snr_lin);
        p_diff_xor(2*ii) = p0_xor_Q-p1_xor_Q;
        logp_xor(2*ii) = log(p0_xor_Q-p1_xor_Q);
        pX_0(2*ii) = p0_xor_Q;
        pX_1(2*ii) = p1_xor_Q;
               
    end
    
    ABits = -1*reshape([pA_0,pA_1]',1,length(pA_0)*2);
    BBits = -1*reshape([pB_0,pB_1]',1,length(pB_0)*2);
    XBits = -1*reshape([pX_0,pX_1]',1,length(pX_0)*2);
end
end

