function [ tx_a, tx_b ] = phy_mod(codeword_a, codeword_b, modulation)

if strcmp(modulation,'BPSK') == 1
    tx_a = 2*codeword_a-1;  %0 map to -1; 1 map to 1
    tx_b = 2*codeword_b-1;
    
elseif strcmp(modulation,'QPSK') == 1
    I_a = codeword_a(1:2:length(codeword_a)-1);
    Q_a = codeword_a(2:2:length(codeword_a));
    tx_a = (2*I_a-1)+1j*(2*Q_a-1);
    tx_a = tx_a/sqrt(2);
    
    I_b = codeword_b(1:2:length(codeword_b)-1);
    Q_b = codeword_b(2:2:length(codeword_b));
    tx_b = (2*I_b-1)+1j*(2*Q_b-1);
    tx_b = tx_b/sqrt(2);
end

end

