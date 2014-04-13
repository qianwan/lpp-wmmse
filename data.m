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

epsilon = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
sM128 = [52, 137, 367, 878, 1819, 3517, 5669];
pM128b24 = [3.2, 15.3, 57.1, 184, 593, 1478, 4500];
pM128b12 = [];
pM128b8 = [2.9, 10.1, 36, 106, 289, 764, 1441];
pM128b6 = [2.7, 9.1, 29.5, ];

sM64 = [51, 134, 373, 908, 1877, 3719, 5120];
pM64b12 = [3, 13.6, 50.8, 168, 434, 897, 1576];
pM64b6 = [];
pM64b5 = [2.9, 10.2, 35.3, ];
pM64b4 = [2.8, 9.6, 32, 90, 214, 425, 689];
semilogx(epsilon, sM128, '-o', epsilon, pM128b24, '-s', epsilon, pM128b8, '-p', ...
         epsilon, sM64, '-^', epsilon, pM64b12, '-d', epsilon, pM64b4, '-h');
legend('Serial BCD, M=128', 'Parallel BCD, M=128, \beta_{ik}=24', 'Parallel BCD, M=128, \beta_{ik}=8', ...
       'Serial BCD, M=64', 'Parallel BCD, M=64, \beta_{ik}=12', 'Parallel BCD, M=64, \beta_{ik}=4');
xlabel('Stopping Criteria');
ylabel('Number of Itertations');
grid on;
