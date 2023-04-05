function [ okX] = PNC_phy_sim(SNR_A, SNR_B, aa, bb)
%matlabpool local 12
%trellis information
constraint_length = 7;

trellis_encoded = poly2trellis(constraint_length,[133,171]);
trellis_decoded = poly2trellis(constraint_length,[133,171]);

modulation = 'BPSK';
%modulation = 'QPSK';

snr = SNR_A;

pkt_length = 64*48/2;  
num_packet = 1;
%degreeB = degree_all(de);
% berA_all = zeros(num_packet,1);
% berB_all = zeros(num_packet,1);
berX_all = zeros(num_packet,1);
% phy_raw_map = zeros(num_packet,3); % collect PHY-layer statistics
for num = 1:num_packet
    %fprintf('packet = %d\n',num);
    degreeB = randi([0,359]);
    msg_a = randi([0 1],pkt_length-8,1); % message sequence for A
    msg_a = [msg_a;zeros(8,1)];
    msg_b = randi([0 1],pkt_length-8,1); % message sequence for B
    msg_b = [msg_b;zeros(8,1)];
    msg_xor = bitxor(msg_a,msg_b);

    codeword_a = convenc(msg_a,trellis_encoded);
    codeword_b = convenc(msg_b,trellis_encoded);

    %modulation
    [ tx_a, tx_b ] = phy_mod(codeword_a, codeword_b, modulation);

    %channel
    n = 1/sqrt(2)*(randn(length(tx_a),1) + j*randn(length(tx_a),1)); 
    
    ha=ones(size(tx_a))*aa;
    tx_ha = ha.*tx_a;

    hbb=ones(size(tx_b)).*bb;
    hb=hbb*exp(j*pi*degreeB/180);
    tx_hb = hb.*tx_b;

%     ha=ones(size(tx_a));
%     tx_ha = ha.*tx_a;
% 
%     hb=ones(size(tx_b))*exp(j*pi*degreeB/180);
%     tx_hb = hb.*tx_b;

    % receiver starts
    rx = tx_ha + sqrt(10^(-snr/10 + SNR_B/10)).*tx_hb + 10^(-snr/20)*n; 

    %demodulation
    [ logp_A, logp_B, logp_xor] = NCMA_demod(rx, ha, hb, modulation, pkt_length, snr);

%     decoded_source_A = viterbi_decoder_bit_level(logp_A,trellis_decoded,pkt_length);
%     decoded_source_B = viterbi_decoder_bit_level(logp_B,trellis_decoded,pkt_length);
    decoded_source_xor = viterbi_decoder_bit_level(logp_xor,trellis_decoded,pkt_length);

%     berA = sum(bitxor(decoded_source_A,msg_a))/length(msg_a);
%     berB = sum(bitxor(decoded_source_B,msg_b))/length(msg_b);
    berX = sum(bitxor(decoded_source_xor,msg_xor))/length(msg_xor);

%     berA_all(num) = berA;
%     berB_all(num) = berB;
    berX_all(num) = berX;

%     okA=0; okB=0; 
    okX=0;
%     if berA == 0
%         okA=1;
%     end
%     if berB == 0
%         okB=1;
%     end
    if berX == 0
        okX = 1;
    end
%     phy_raw_map(num,:) = [okA, okB, okX];

end




%matlabpool close
