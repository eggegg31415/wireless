%% BPSK transmission over AWGN channel
close all;clear all;clc;           
dist=100:100:400;       % distance in meters
PtdBm=10;               % transmit power in dBm
PndBm=-85;              % noise power in dBm
Pt=10^(PtdBm/10)/1000;  % transmit power in watt
Pn=10^(PndBm/10)/1000;  % noise power in watt
Bit_Length=1e3;         % number of bits transmitted
MODORDER = [1,2,4];     % modulation orders

%% Friss Path Loss Model
Gt=1;
Gr=1;
freq=2.4e9;
lambda=3e8/freq;
Pr=Pt*Gt*Gr*(lambda./(4*pi*dist)).^2;
PrdBm=log10(Pr*1000)*10;
SNRdB=PrdBm - PndBm;
SNR=10.^(SNRdB/10);
NumStream = 2;  % MIMO: Number of streams

%% Generate bit streams
tx_data = randi(2, 1, Bit_Length) - 1;          

% MIMO: update NumSym
NumSym(MODORDER) = length(tx_data)./MODORDER;

%% Constellation points
% BPSK: {1,0} -> {1+0i, -1+0i}
% QPSK: {11,10,01,00} -> {1+i, -1+i, -1-i, 1-i} * scaling factor
% 16QAM: {1111,1110,1101,1100,1011,1010,1001,1000,0111,0110,0101,0100,0011,0110,0001,0000}
% -> {3a+3ai,3a+ai,a+3ai,a+ai,-a+3ai,-3a+3ai,-3a+ai,3a-ai,3a-3ai,a-ai,a-3i,-a-ai,-a-3ai,-3a-ai,-3a-3ai}

BPSKBit = [0; 1];
BPSK = [-1+0i; 1+0i];
QPSKBit = [0 0; 0 1; 1 0; 1 1];
QPSK = [1-i, -1-i, -1+i, 1+i]./sqrt(2);
QAMBit = [1 1 1 1; 1 1 1 0; 1 1 0 1; 1 1 0 0; 1 0 1 1; 1 0 1 0; 1 0 0 1; 1 0 0 0; 0 1 1 1; 0 1 1 0; 0 1 0 1; 0 1 0 0; 0 0 1 1; 0 0 1 0; 0 0 0 1; 0 0 0 0];
QAM = [3+3i, 3+i, 1+3i, 1+1i, -1+3i, -1+i, -3+3i, -3+i, 3-i, 3-3i, 1-i, 1-3i, -1-i, -1-3i, -3-i, -3-3i]./sqrt(10);
IQPoint(4,:) = QAM;
IQPoint(2,1:4) = QPSK;
IQPoint(1,1:2) = BPSK;

n=(randn(NumStream,Bit_Length)+randn(NumStream, Bit_Length)*i)/sqrt(2);  % MIMO: AWGN noises
n=n*sqrt(Pn);

% repeat 5 times
for round = 1:5
    %% MIMO channel: h dimension:  NumStream x NumStream
    h = (randn(NumStream, NumStream) + randn(NumStream, NumStream) * 1i);
    h = h ./ abs(h);
    
    % update theta
    correlation = abs(real(dot(h(:, 1).', h(:, 2)))) / (norm(h(:, 1)) * norm(h(:, 2)));
    theta(round) = rad2deg(acos(correlation));
    
    w = inv(h);
    amp(1, round) = norm(abs(w(1, :)))^2;
    amp(2, round) = norm(abs(w(2, :)))^2;
    
    for mod_order = MODORDER

        %% modulation
        if (mod_order == 1)
            % BPSK
            [ans ix] = ismember(tx_data', BPSKBit, 'rows'); 
            s = BPSK(ix).';
        elseif (mod_order == 2)
            % QPSK
            tx_data_reshape = reshape(tx_data, length(tx_data)/mod_order, mod_order);
            [ans ix] = ismember(tx_data_reshape, QPSKBit, 'rows');
            s = QPSK(ix);
        else
            % QAM
            tx_data_reshape = reshape(tx_data, length(tx_data)/mod_order, mod_order);
            [ans ix] = ismember(tx_data_reshape, QAMBit, 'rows');
            s = QAM(ix);
        end

        % MIMO: reshape to NumStream streams
        x = reshape(s, NumStream, length(s)/NumStream);

        for d=1:length(dist)
            
            %% transmission with noise
            % update Y = HX + N
            y = h * x * sqrt(Pr(d)) + n(:,1:(Bit_Length/(2*mod_order)));

            %% ZF equalization
            x_est = w * y / sqrt(Pr(d));
            s_est = reshape(x_est, [1,(Bit_Length/mod_order)]);
            re_x = reshape(x, [1,Bit_Length/mod_order]);

            %% demodulation
            tx_predict = demodulation(mod_order, s_est, Bit_Length);
            tx_predict = reshape(tx_predict, mod_order, (Bit_Length/mod_order));
            tx_predict = tx_predict';
            tx_predict = reshape(tx_predict, 1, Bit_Length);

            error(d,round,mod_order) = 0;
            for i=1:Bit_Length
                if tx_predict(i) ~= tx_data(i) error(d,round,mod_order) = error(d,round,mod_order)+1;, end;
            end
            
            SNR(round,d,mod_order)=Pr(d)/Pn;
            SNRdB(round,d,mod_order)=10*log10(SNR(round,d,mod_order));
            
            BER_simulated(d,round,mod_order) = error(d,round,mod_order) / Bit_Length;
            noise = s_est - re_x;
            SNR(round,d,mod_order)= 1 / mean( abs(noise.^2) );
            SNRdB_simulated(d,round,mod_order)=10*log10(SNR(round,d,mod_order));

            % Uncomment to show the plot of y
            %{
            subplot(2, 2, d)
            hold on;

            plot(s_est,'bx'); 
            plot(s,'ro');
            hold off;
            xlim([-2,2]);
            ylim([-2,2]);
            title(sprintf('Constellation points d=%d', dist(d)));
            legend('decoded samples', 'transmitted samples');
            grid
            %}
        end
        % filename = sprintf('IQ_%d.jpg', mod_order);
        % saveas(gcf,filename,'jpg')
    end
end

figure('units','normalized','outerposition',[0 0 1 1])
hold on;
bar(dist,SNRdB_simulated(:,:,1));
plot(dist,SNRdB(1,:,1),'bx-', 'Linewidth', 1.5);
hold off;
title('SNR');
xlabel('Distance [m]');
ylabel('SNR [dB]');
legend('simu-1', 'simu-2', 'simu-3', 'simu-4', 'simu-5', 'siso-theory');
axis tight 
grid
saveas(gcf,'SNR.jpg','jpg')

figure('units','normalized','outerposition',[0 0 1 1])
hold on;
bar(1:5, theta);
hold off;
title('channel angle');
xlabel('Iteration index');
ylabel('angle [degree]');
axis tight 
grid
saveas(gcf,'angle.jpg','jpg')

figure('units','normalized','outerposition',[0 0 1 1])
hold on;
bar(1:5, amp);
hold off;
title('Amplification');
xlabel('Iteration index');
ylabel('noise amplification');
legend('x1', 'x2');
axis tight 
grid
saveas(gcf,'amp.jpg','jpg')

figure('units','normalized','outerposition',[0 0 1 1])
hold on;
plot(dist,mean(BER_simulated(:,:,1)',1),'bo-','linewidth',2.0);
plot(dist,mean(BER_simulated(:,:,2)',1),'rv--','linewidth',2.0);
plot(dist,mean(BER_simulated(:,:,4)',1),'mx-.','linewidth',2.0);
hold off;
title('BER');
xlabel('Distance [m]');
ylabel('BER');
legend('BPSK','QPSK','16QAM');
axis tight 
grid
saveas(gcf,'BER.jpg','jpg')
return;

function demod = demodulation(mod_order, s, Bit_Length)
    len = Bit_Length/mod_order;
    demod = [];
    if mod_order == 1
        for i = 1 : len
            if real(s(i)) >= 0 demod(i) = 1;, else demod(i) = 0;, end;
        end
    elseif mod_order == 2  
        scale = sqrt(2);
        s = s * scale;
        QPSK = [1 -1; -1 -1; -1 1; 1 1];
        QPSKBit = [0 0; 0 1; 1 0; 1 1];        
        T = [1 2 3; 4 1 2;];
        TR = triangulation(T, QPSK);
        for i = 1 : len            
            [ID,d] = nearestNeighbor(TR, [real(s(i)) imag(s(i))]);
            demod = [demod QPSKBit(ID, :)];
        end
    elseif mod_order == 4
        scale = sqrt(10);
        s = s * scale;
        QAM = [3 3; 3 1; 1 3; 1 1;
            -1 3; -1 1; -3 3; -3 1; 
            3 -1; 3 -3; 1 -1; 1 -3; 
            -1 -1; -1 -3; -3 -1; -3 -3];
        QAMBit = [1 1 1 1; 1 1 1 0; 1 1 0 1; 1 1 0 0; 
                1 0 1 1; 1 0 1 0; 1 0 0 1; 1 0 0 0; 
                0 1 1 1; 0 1 1 0; 0 1 0 1; 0 1 0 0; 
                0 0 1 1; 0 0 1 0; 0 0 0 1; 0 0 0 0];    
        T = [1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 1 2;];
        TR = triangulation(T, QAM);
        for i = 1 : len
            [ID,d] = nearestNeighbor(TR, [real(s(i)) imag(s(i))]);
            demod = [demod QAMBit(ID, :)];
        end
    end
end
