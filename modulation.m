%% BPSK transmission over AWGN channel
close all;clear all;clc;    % BPSK
dist=100:100:400;           % distance in meters
PtdBm=10;                   % transmit power in dBm
PndBm=-85;                  % noise power in dBm
Pt=10^(PtdBm/10)/1000;      % transmit power in watt
Pn=10^(PndBm/10)/1000;      % noise power in watt
Bit_Length=1e3;             % number of bits transmitted

%% Friss Path Loss Model
Gt=1;
Gr=1;
freq=2.4e9;

light = 3e8;
lambda = light/freq;

Pr = (Pt * Gt * Gr * ((lambda./(4*pi*dist)).^2))';

%% AWGN channel
tx_data = randi(2, 1, Bit_Length) - 1;                  % random between 0 and 1
n=(randn(1,Bit_Length)+randn(1,Bit_Length)*i)/sqrt(2);  % AWGN noises
n=n*sqrt(Pn);

mod_order_list = [1,2,4];
for mod_order=mod_order_list
    x(mod_order, :) = modulation(mod_order, tx_data);

    for d=1:length(dist)
        y(mod_order,d,:)=sqrt(Pr(d))*x(mod_order,:)+n;
    end
end

%% Equalization
% Detection Scheme:(Soft Detection)
% Error if input and output are of different signs
for mod_order=mod_order_list
    figure('units','normalized','outerposition',[0 0 1 1])
    sgtitle(sprintf('Modulation order: %d', mod_order));
    len = Bit_Length/mod_order;
    for d=1:length(dist)
        s = reshape((y(mod_order, d, :)./sqrt(Pr(d))), [1, Bit_Length]);
        err = 0;
        tx_predict = demodulation(mod_order, s);
        for i=1:Bit_Length
            if tx_predict(i) ~= tx_data(i) err = err+1;, end;
        end
        noise = s - x(mod_order, :);

        BER_simulated(mod_order, d)=err/Bit_Length;
        SNR_simulated(mod_order, d) = 1/mean(abs(noise(1:len).^2));
        SNRdB_simulated(mod_order, d) = 10*log10(SNR_simulated(mod_order, d));

        subplot(2, 2, d)
        hold on;
        plot(s(1:length(s)/mod_order),'bx');
        if mod_order == 1
            plot(x(mod_order, :),0i,'ro');
        else
            plot(x(mod_order, 1:length(x)/mod_order),'ro');
        end
        hold off;
        xlim([-2,2]);
        ylim([-2,2]);
        title(sprintf('Constellation points d=%d', dist(d)));
        legend('decoded samples', 'transmitted samples');
        grid
    end
    filename = sprintf('IQ_%d.jpg', mod_order);
    saveas(gcf,filename,'jpg')
end

figure('units','normalized','outerposition',[0 0 1 1])
hold on;
semilogy(dist,SNRdB_simulated(1, :),'bo-','linewidth',2.0);
semilogy(dist,SNRdB_simulated(2, :),'rv--','linewidth',2.0);
semilogy(dist,SNRdB_simulated(4, :),'mx-.','linewidth',2.0);
hold off;
title('SNR');
xlabel('Distance [m]');
ylabel('SNR [dB]');
legend('BPSK','QPSK','16QAM');
axis tight
grid
saveas(gcf,'SNR.jpg','jpg')

figure('units','normalized','outerposition',[0 0 1 1])
hold on;
semilogy(dist,BER_simulated(1, :),'bo-','linewidth',2.0);
semilogy(dist,BER_simulated(2, :),'rv--','linewidth',2.0);
semilogy(dist,BER_simulated(4, :),'mx-.','linewidth',2.0);
hold off;
title('BER');
xlabel('Distance [m]');
ylabel('BER');
legend('BPSK','QPSK','16QAM');
axis tight
grid
saveas(gcf,'BER.jpg','jpg')
return;

function mod = modulation(mod_order, data)
    len = length(data)/mod_order;
    mod = zeros(1, length(data));
    if mod_order == 1
        %% BPSK: {1,0} -> {1+0i, -1+0i}
        for i = 1 : len
            mod(i) = 2*data(i) - 1 ;
        end
    elseif mod_order == 2
        %% QPSK: {11,10,01,00} -> {1+i, -1+i, -1-i, 1-i} * scaling factor
        cnt = 1;
        scale = 1/sqrt(2);
        QPSKBit = [1 1; 1 0; 0 1; 0 0];
        QPSK = [1+j, -1+j, -1-j, 1-j];
        for i = 1 : len
            pos = ismember(QPSKBit, data(cnt:cnt+1), "rows");
            mod(i) = QPSK(find(pos));
            cnt = cnt + 2;
        end
        mod = mod*scale;
    elseif mod_order == 4
        %% 16QAM: {1111, 1110, 1101, 1100, 1011, 1010, 1001, 1000, 0111, 0110, 0101, 0100, 0011, 0010, 0001, 0000}
        %% -> {3a+3ai, 3a+ai, a+3ai, a+ai, -a+3ai, -3a+3ai, -3a+ai, -a+ai, 3a-ai, 3a-3ai, a-ai, a-3i, -a-ai, -a-3ai, -3a-ai, -3a-3ai}
        cnt=1;
        scale = sqrt(4/39);
        QAMBit = [1 1 1 1; 1 1 1 0; 1 1 0 1; 1 1 0 0; 
                1 0 1 1; 1 0 1 0; 1 0 0 1; 1 0 0 0; 
                0 1 1 1; 0 1 1 0; 0 1 0 1; 0 1 0 0; 
                0 0 1 1; 0 0 1 0; 0 0 0 1; 0 0 0 0];
        QAM = [3+3j, 3+j, 1+3j, 1+1j, ...,
            -1+3j, -3+3j, -3+j, -1+j, ...,
            3-j, 3-3j, 1-j, 1-3j, ...,
            -1-j, -1-3j, -3-j, -3-3j];
        for i = 1 : len
            pos = ismember(QAMBit, data(cnt:cnt+3), "rows");
            mod(i) = QAM(find(pos));
            cnt = cnt + 4;
        end
        mod = mod*scale;
    end
end

function demod = demodulation(mod_order, s)
    len = length(s)/mod_order;
    demod = zeros(1, length(s));
    if mod_order == 1
        %% BPSK: {1,0} -> {1+0i, -1+0i}
        for i = 1 : len
            if real(s(i)) >= 0 demod(i) = 1;, else demod(i) = 0;, end;
        end
    elseif mod_order == 2
        %% QPSK: {11,10,01,00} -> {1+i, -1+i, -1-i, 1-i} * scaling factor
        cnt = 1;
        for i = 1 : len
            if real(s(i)) >= 0 && imag(s(i)) >= 0
                demod(cnt) = 1; demod(cnt+1) = 1;
            elseif real(s(i)) < 0 && imag(s(i)) >= 0
                demod(cnt) = 1; demod(cnt+1) = 0;
            elseif real(s(i)) < 0 && imag(s(i)) < 0
                demod(cnt) = 0; demod(cnt+1) = 1;
            elseif real(s(i)) >= 0 && imag(s(i)) < 0
                demod(cnt) = 0; demod(cnt+1) = 0;
            end
            cnt = cnt + 2;
        end
    elseif mod_order == 4
        %% 16QAM: {1111, 1110, 1101, 1100, 1011, 1010, 1001, 1000, 0111, 0110, 0101, 0100, 0011, 0010, 0001, 0000}
        %% -> {3a+3ai, 3a+ai, a+3ai, a+ai, -a+3ai, -3a+3ai, -3a+ai, -a+ai, 3a-ai, 3a-3ai, a-ai, a-3i, -a-ai, -a-3ai, -3a-ai, -3a-3ai}
        cnt=1;
        scale = sqrt(4/39);
        s = s/scale;
        for i = 1 : len
            if real(s(i)) >= 2 && imag(s(i)) > 2
                demod(cnt) = 1; demod(cnt+1) = 1; demod(cnt+2) = 1; demod(cnt+3) = 1;
            elseif real(s(i)) >= 2 && imag(s(i)) > 0 && imag(s(i)) < 2
                demod(cnt) = 1; demod(cnt+1) = 1; demod(cnt+2) = 1; demod(cnt+3) = 0;
            elseif real(s(i)) >= 0 && real (s(i)) < 2 && imag(s(i)) > 2
                demod(cnt) = 1; demod(cnt+1) = 1; demod(cnt+2) = 0; demod(cnt+3) = 1;
            elseif real(s(i)) >= 0 && real (s(i)) < 2 && imag(s(i)) > 0 && imag(s(i)) < 2
                demod(cnt) = 1; demod(cnt+1) = 1; demod(cnt+2) = 0; demod(cnt+3) = 0;
            elseif real(s(i)) >= -2 && real(s(i)) < 0 && imag(s(i)) > 2
                demod(cnt) = 1; demod(cnt+1) = 0; demod(cnt+2) = 1; demod(cnt+3) = 1;
            elseif real(s(i)) < -2 && imag(s(i)) > 2
                demod(cnt) = 1; demod(cnt+1) = 0; demod(cnt+2) = 1; demod(cnt+3) = 0;
            elseif real(s(i)) < -2 && imag(s(i)) > 0 && imag(s(i)) < 2
                demod(cnt) = 1; demod(cnt+1) = 0; demod(cnt+2) = 0; demod(cnt+3) = 1;
            elseif real(s(i)) >= -2 && real(s(i)) < 0 && imag(s(i)) > 0 && imag(s(i)) < 2
                demod(cnt) = 1; demod(cnt+1) = 0; demod(cnt+2) = 0; demod(cnt+3) = 0;
            elseif real(s(i)) >= 2 && imag(s(i)) < 0 && imag(s(i)) > -2
                demod(cnt) = 0; demod(cnt+1) = 1; demod(cnt+2) = 1; demod(cnt+3) = 1;
            elseif real(s(i)) >= 2 && imag(s(i)) < -2
                demod(cnt) = 0; demod(cnt+1) = 1; demod(cnt+2) = 1; demod(cnt+3) = 0;
            elseif real(s(i)) >= 0 && real (s(i)) < 2 && imag(s(i)) < 0 && imag(s(i)) > -2
                demod(cnt) = 0; demod(cnt+1) = 1; demod(cnt+2) = 0; demod(cnt+3) = 1;
            elseif real(s(i)) >= 0 && real (s(i)) < 2 && imag(s(i)) < -2
                demod(cnt) = 0; demod(cnt+1) = 1; demod(cnt+2) = 0; demod(cnt+3) = 0;
            elseif real(s(i)) >= -2 && real(s(i)) < 0 && imag(s(i)) < 0 && imag(s(i)) > -2
                demod(cnt) = 0; demod(cnt+1) = 0; demod(cnt+2) = 1; demod(cnt+3) = 1;
            elseif real(s(i)) >= -2 && real(s(i)) < 0 && imag(s(i)) < -2
                demod(cnt) = 0; demod(cnt+1) = 0; demod(cnt+2) = 1; demod(cnt+3) = 0;
            elseif real(s(i)) < -2 && imag(s(i)) < 0 && imag(s(i)) > -2
                demod(cnt) = 0; demod(cnt+1) = 0; demod(cnt+2) = 0; demod(cnt+3) = 1;
            elseif real(s(i)) < -2 && imag(s(i)) < -2
                demod(cnt) = 0; demod(cnt+1) = 0; demod(cnt+2) = 0; demod(cnt+3) = 0;
            end
            cnt = cnt + 4;
        end
    end
end

