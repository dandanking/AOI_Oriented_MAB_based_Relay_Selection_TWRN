function [ ok] = single_phy_sim(SNR, aa)
constraint_length = 7;

trellis_encoded = poly2trellis(constraint_length,[133,171]);
trellis_decoded = poly2trellis(constraint_length,[133,171]);

modulation = 'BPSK';
snr = SNR;

pkt_length = 64*48/2;
num_packet = 1;
berA_all = zeros(num_packet,1);
phy_raw_map = zeros(num_packet,1); % collect PHY-layer statistics

for num = 1:num_packet
%         fprintf('packet = %d\n',num);
    msg_a = randi([0 1],pkt_length-8,1); % message sequence for A
    msg_a = [msg_a;zeros(8,1)];

    codeword_a = convenc(msg_a,trellis_encoded);

    %modulation
    [ tx_a] = phy_mod(codeword_a, 0,modulation);

    %channel
    n = 1/sqrt(2)*(randn(length(tx_a),1) + 1i*randn(length(tx_a),1)); 

    ha=ones(size(tx_a))*aa;
    tx_ha = ha.*tx_a;               

    % receiver starts
    rx = tx_ha + 10^(-snr/20)*n; 

    %demodulation
    [ logp_A,~] = single_demod(rx, ha,modulation, pkt_length, snr);
    decoded_source_A = viterbi_decoder_bit_level(logp_A,trellis_decoded,pkt_length);
    berA = sum(bitxor(decoded_source_A,msg_a))/length(msg_a);
    berA_all(num) = berA;
    ok=0;
    if berA == 0
        ok=1;
    end
    phy_raw_map(num,:) = ok;
end



%matlabpool close
