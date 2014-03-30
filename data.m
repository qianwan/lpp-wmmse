dB = [0 : 5 : 30];
sw = [200 300 430 580 680 760 765];
pw = [227 356 496 607 724 767 802];
ww = [280 390 510 680 790 885 910];
plot(dB, sw, '^-', dB, pw, 'o-', dB, ww, 's-');
legend('S-WMMSE', 'LPS-WMMSE', 'WMMSE');
xlabel('SNR (dB)');
ylabel('System Throughput (bits per channel use)')
grid on;

dB = [0 : 5 : 30];
ssv = [1 2.4 6 8.2 11.8 12.5 12.45];
psv = [0.69 1.8 5.1 6.7 10.1 11.3 11.7];
wsv = [6.1 13.1 18 19.1 19.5 19.5 19.5];
plot(dB, ssv, '^-', dB, psv, 'o-', dB, wsv, 's-');
legend('S-WMMSE', 'LPS-WMMSE', 'WMMSE');
xlabel('SNR (dB)');
ylabel('Average Number of Serving BSs per User');
grid on;

epsilon = [1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7];
iterSer = [0.39 8.4 24.6 58.4 122 202.5 308.8];
iterPal = [0.69 4.5 12.3 27.6 48.1 78.3 118.2];
iterPal4 = [0.65 3.9 9.9 21.0 36.1 55.4 82.7];
semilogx(epsilon, iterSer, '^-', epsilon, iterPal, 'o-', epsilon, iterPal4, 's-');
legend('Serial BCD', 'Parallel BCD, \beta=6', 'Parallel BCD, \beta=4');
xlabel('Stopping criteria');
ylabel('Number of iterations');
grid on;
