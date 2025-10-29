# dc-lab1
all codes
%EXPT-1
clc; clear;
V = [1 0 0; 1 1 1; 0 0 1]';
num_vectors = size(V, 2);
U = zeros(size(V));     
E = zeros(size(V));     
U(:,1) = V(:,1);
E(:,1) = U(:,1) / norm(U(:,1));
for i = 2:num_vectors
    U(:,i) = V(:,i);
    for j = 1:i-1
        U(:,i) = U(:,i) - (dot(E(:,j), V(:,i))) * E(:,j);
    end
    E(:,i) = U(:,i) / norm(U(:,i));
end
disp('Orthonormal basis vectors (columns):');
disp(E);
figure; hold on; grid on; axis equal;
quiver3(0,0,0, V(1,1), V(2,1), V(3,1), 'r','LineWidth',2);
quiver3(0,0,0, V(1,2), V(2,2), V(3,2), 'g','LineWidth',2);
quiver3(0,0,0, V(1,3), V(2,3), V(3,3), 'b','LineWidth',2);
quiver3(0,0,0, E(1,1), E(2,1), E(3,1), 'r--','LineWidth',2);
quiver3(0,0,0, E(1,2), E(2,2), E(3,2), 'g--','LineWidth',2);
quiver3(0,0,0, E(1,3), E(2,3), E(3,3), 'b--','LineWidth',2);
xlabel('X'); ylabel('Y'); zlabel('Z');
legend('v1','v2','v3','e1','e2','e3');
title('Original and Orthonormal Basis Vectors');

expt-2
clc;
clear all;
close all;
bit_seq = [1 1 0 0 0 0 1 1];
N = length(bit_seq);
fc = 1;
t = 0:0.001:2; 

b = [];    
qpsk1 = []; 
bec = [];   
bes = [];  
bit_e = []; 
bit_o = []; 


for i = 1:N
    bx = bit_seq(i) * ones(1, 1000);
    b = [b, bx];
end


bit_seq(bit_seq == 0) = -1;


b_o = bit_seq(1:2:end); 
b_e = bit_seq(2:2:end); 


for i = 1:length(b_e)
    be_c = b_e(i) * cos(2*pi*fc*t); % Even bits on cosine
    bo_s = b_o(i) * sin(2*pi*fc*t); % Odd bits on sine
    q = be_c + bo_s;                % Combined QPSK signal
    
    qpsk1 = [qpsk1, q];
    bec = [bec, be_c];
    bes = [bes, bo_s];
    
    bit_e = [bit_e, b_e(i)*ones(1, length(t))];
    bit_o = [bit_o, b_o(i)*ones(1, length(t))];
end


figure('Name', 'QPSK Modulation', 'NumberTitle', 'off');

subplot(5,1,1);
plot(b, 'LineWidth', 1.5);
grid on;
axis([0 N*1000 -0.5 1.5]);
title('Binary Input Sequence');
ylabel('Amplitude');

subplot(5,1,2);
plot(bes, 'b', 'LineWidth', 1.2); hold on;
plot(bit_o, 'r--');
grid on;
axis([0 N*1000 -1.5 1.5]);
title('Odd Bits (Sine Component)');
ylabel('Amplitude');

subplot(5,1,3);
plot(bec, 'g', 'LineWidth', 1.2); hold on;
plot(bit_e, 'r--');
grid on;
axis([0 N*1000 -1.5 1.5]);
title('Even Bits (Cosine Component)');
ylabel('Amplitude');

subplot(5,1,4);
plot(qpsk1, 'k', 'LineWidth', 1.2);
grid on;
axis([0 N*1000 -2 2]);
title('QPSK Modulated Signal');
ylabel('Amplitude');


subplot(5,1,5);
constellation = [1+1j, -1+1j, -1-1j, 1-1j];
plot(real(constellation), imag(constellation), 'bo', 'MarkerSize', 8, 'LineWidth', 2);
grid on;
axis([-2 2 -2 2]);
title('QPSK Constellation Diagram');
xlabel('In-phase (I)');
ylabel('Quadrature (Q)');




expt-3
N = 1e4;
SNR_dB = 0:5:20;
pulse_width = 1;
data = randi([0 1], N, 1);
t = 0:0.01:pulse_width;
rect_pulse = ones(size(t));
BER = zeros(length(SNR_dB), 1);
for snr_idx = 1:length(SNR_dB)
tx_signal = [];
for i = 1:N
if data(i) == 1
tx_signal = [tx_signal; rect_pulse'];
else
tx_signal = [tx_signal; zeros(size(rect_pulse'))];
end
end
SNR = 10^(SNR_dB(snr_idx) / 10);
noise_power = 1 / (2 * SNR);
noise = sqrt(noise_power) * randn(length(tx_signal), 1);
rx_signal = tx_signal + noise;
matched_filter = rect_pulse;
filtered_signal = conv(rx_signal, matched_filter, 'same');
sample_interval = round(length(filtered_signal) / N);
sampled_signal = filtered_signal(1:sample_interval:end);
estimated_bits = sampled_signal > 0.5;
num_errors = sum(estimated_bits ~= data);
BER(snr_idx) = num_errors / N;
end
figure;
semilogy(SNR_dB, BER, 'b-o');
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs. SNR for Rectangular Pulse Modulated Bin

expt 4

M= 16;
N= 1000;
bits = randi([01], 1,N);
symbols = zeros(1, N/4);
for i = 1:N/4
symbols(i) = (2*bits(4*i-3)-1) + 1j*(2*bits(4*i-2)-1) + 2*(2*bits(4*i-1)-1) + 2j*(2*bits(4*i)-1);
end

scatter(real(symbols), imag(symbols), 'bo');
grid on;
xlabel('In-phase');
ylabel('Quadrature');
title('16-QAM Constellation');

snr_db=20;
rx_signal=awgn(symbols,snr_db,'measured');
figure;
plot(real(rx_signal),imag(rx_signal),'bo','Markersize',6, 'Linewidth',2)
xlabel('In-phase');
ylabel('Quadrature');
title('16 OAM constellation with noise')
grid on; 
axis equal;
axis([-4 4 -4 4]); 




expt 5
clc; clear;
p = input('Enter the probabilities');
n = length(p);
symbols = 1:n;
[dict, avglen] = huffmandict(symbols, p);
disp('The Huffman code dictionary:');
for i = 1:n
    fprintf('Symbol %d: %s',symbols(i), num2str(dict{i, 2}));
end
sym = input(sprintf('Enter the symbols between 1 to %d in []: ', n));
encoded = huffmanenco(sym, dict);
disp('The encoded output:');
disp(encoded);
bits = input('Enter the bit stream in []: ');
decoded = huffmandeco(bits, dict);
disp('The decoded symbols are:');
disp(decoded);

expt 6

clc;
clear;
data = [1 0 1 0];
p1 = mod(data(1) + data(2) + data(4), 2);
p2 = mod(data(1) + data(3) + data(4), 2);
p3 = mod(data(2) + data(3) + data(4), 2);
encoded_data = [p1 p2 data(1) p3 data(2) data(3) data(4)];
disp('Encoded Data:');
disp(encoded_data);
recieved_data=[1 0 1 1 0 1 0];
disp('recieved data');
disp(recieved_data);
s1 = mod(recieved_data(1) + recieved_data(3) + recieved_data(5) + recieved_data(7), 2);
s2 = mod(recieved_data(2) + recieved_data(3) + recieved_data(6) +recieved_data(7), 2);
s3 = mod(recieved_data(4) + recieved_data(5) + recieved_data(6) + recieved_data(7), 2);
error_location = bin2dec([num2str(s3) num2str(s2) num2str(s1)]);
if error_location ~= 0
    recv(error_location) = mod(recv(error_location) + 1, 2);
    fprintf('Bit error corrected at position %d.', error_location);
else
    disp('No error detected.');
end
decoded_data = recieved_data([3 5 6 7]);
disp('Decoded Data (4 bits):');
disp(decoded_data);







