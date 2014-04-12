dB = [0 : 5 : 30];
rM8Q8C7 = [259, 334, 422, 510, 590, 650, 692];
rpcM8Q8C7 = [1.0 0.99 0.99 0.98 0.99, 0.99, 0.99];
wmmseM8Q8C7 = [278, 357, 430, 495, 536, 562, 584]

rM8Q16C7 = [289, 390, 482, 596, 722, 845, 957];

rM4Q5C7 = [130, 175, 233, 286, 315, 336, 342];
wmmseM4Q5C7 = [164, 213, 258, 290, 309, 315, 319];
% power consumption 100%

plot(dB, rM8Q16C7, '-o', dB, rM8Q8C7, '-s', dB, wmmseM8Q8C7, '-p', dB, rM4Q5C7, '-^', dB, wmmseM4Q5C7, '-*', 'LineWidth', 1.5);
legend('LPP-WMMSE, M=128', 'LPP-WMMSE, M=64', 'WMMSE, M=64', 'LPP-WMMSE, M=20', 'WMMSE, M=20');
grid on;
xlabel('SNR (dB)');
ylabel('System Throughput (bits per channel user)');
