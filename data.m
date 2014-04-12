dB = [0 : 5 : 30];
rM8Q8C7 = [259, 340, 422, 510, 590, 650, 692];
wmmseM8Q8C7 = [278, 357, 430, 495, 536, 562, 584];

rM8Q16C7 = [289, 390, 482, 596, 722, 845, 957];

rM4Q5C7 = [130, 175, 233, 286, 315, 336, 342];
wmmseM4Q5C7 = [164, 213, 258, 290, 309, 315, 319];
% power consumption 100%

plot(dB, rM8Q16C7, '-o', dB, rM8Q8C7, '-s', dB, wmmseM8Q8C7, '-p', dB, rM4Q5C7, '-^', dB, wmmseM4Q5C7, '-d', 'LineWidth', 1.5);
legend('LPP-WMMSE, M=128, G=16', 'LPP-WMMSE, M=64, G=8', 'WMMSE, M=64', 'LPP-WMMSE, M=20, G=5', 'WMMSE, M=20');
grid on;
xlabel('SNR (dB)');
ylabel('System Throughput (bits per channel user)');

dB = [5 : 5 : 30];
rpcWMMSE = [1, 1, 1, 1, 1, 1];
rpcM8Q8C7 = [0.99 0.99 0.98 0.96, 0.93, 0.905];
rpcM8Q16C7 = [0.95 0.90 0.87 0.855 0.83 0.81];
bar1 = bar(dB, rpcM8Q8C7, 'BarWidth', 0.2, 'FaceColor', 'b');
hold on;
bar2 = bar(dB - 1, rpcWMMSE, 'BarWidth', 0.2, 'FaceColor', 'k');
hold on;
bar3 = bar(dB + 1, rpcM8Q16C7, 'BarWidth', 0.2, 'FaceColor', 'r');
xlabel('SNR (dB)');
ylabel('Relative Power Consumption');
grid on;
legend('LPP-WMMSE, M=64', 'WMMSE', 'LPP-WMMSE, M=128');