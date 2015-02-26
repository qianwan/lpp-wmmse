dB = [0 : 5 : 30];
rM8Q8C7 = [259, 340, 422, 510, 590, 650, 692];
wmmseM8Q8C7 = [278, 357, 430, 495, 536, 562, 584];

rM8Q16C7 = [289, 390, 482, 595, 719, 837, 957];
wmmseM8Q16C7 = [300, 370, 445, 540, 635, 720, 825];

rM4Q5C7 = [130, 175, 233, 286, 315, 336, 342];
wmmseM4Q5C7 = [164, 213, 258, 290, 309, 315, 319];
% power consumption 100%

plot(dB, rM8Q16C7, '-o', dB, wmmseM8Q16C7, '-s', dB, rM8Q8C7, '-p', dB, wmmseM8Q8C7, '-^', dB, rM4Q5C7, '-d', dB, wmmseM4Q5C7, '-h', 'LineWidth', 1.5);
legend('LPP-WMMSE, M=128, G=16', 'WMMSE, M=128', 'LPP-WMMSE, M=64, G=8', 'WMMSE, M=64', 'LPP-WMMSE, M=20, G=5', 'WMMSE, M=20');
grid on;
xlabel('SNR (dB)');
ylabel('系统吞吐率(比特每用户)');

dB = [5 : 5 : 30];
rpcWMMSE = [1, 1, 1, 1, 1, 1];
rpcM8Q8C7 = [0.99 0.99 0.98 0.96, 0.93, 0.905];
rpcM8Q16C7 = [0.95 0.90 0.87 0.855 0.83 0.81];
bar1 = bar(dB, rpcM8Q8C7, 'BarWidth', 0.2, 'FaceColor', 'b');
hold on;
bar2 = bar(dB - 1, rpcWMMSE, 'BarWidth', 0.2, 'FaceColor', 'c');
hold on;
bar3 = bar(dB + 1, rpcM8Q16C7, 'BarWidth', 0.2, 'FaceColor', 'r');
xlabel('SNR (dB)');
ylabel('功率消耗相对值');
grid on;
legend('LPP-WMMSE, M=64', 'WMMSE', 'LPP-WMMSE, M=128');
axis([0 35 0.6 1]);

range = [1:7];
epsilon = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
%epsilon = [1, 2, 3, 4, 5, 6, 7];
sM128 = [52, 137, 367, 878, 1819, 3517, 5669];
pM128b24 = [3.2, 15.3, 57.1, 184, 593, 1478, 3189];
pM128b12 = [];
pM128b8 = [2.9, 10.1, 36, 106, 289, 764, 1441];
pM128b6 = [2.7, 9.1, 29.5, ];

sM64 = [51, 134, 373, 908, 1877, 3719, 5120];
pM64b12 = [3, 13.6, 50.8, 168, 434, 897, 1576];
pM64b6 = [];
pM64b5 = [2.9, 10.2, 35.3, ];
pM64b4 = [2.8, 9.6, 32, 90, 214, 425, 689];
loglog(epsilon(range), sM128(range), '-o', epsilon(range), pM128b24(range), '-s', epsilon(range), pM128b8(range), '-p', ...
         epsilon(range), sM64(range), '-^', epsilon(range), pM64b12(range), '-d', epsilon(range), pM64b4(range), '-h', 'LineWidth', 1.5);
legend('序列BCD, M=128', '并行BCD, M=128, \beta_{ik}=24', '并行BCD, M=128, \beta_{ik}=8', ...
       '序列BCD, M=64', '并行BCD, M=64, \beta_{ik}=12', '并行BCD, M=64, \beta_{ik}=4');
xlabel('停止条件');
ylabel('迭代次数');
grid on;



dB = [0 : 5 : 30];
SWMMSE   = [200, 300, 429, 580, 682, 760, 770];
LPPWMMSE = [228, 355, 495, 605, 725, 790, 802];
WMMSE    = [282, 392, 504, 677, 787, 885, 904];
plot(dB, SWMMSE, '-^', dB, LPPWMMSE, '-o', dB, WMMSE, '-s', 'LineWidth', 1.5);
legend('S-WMMSE', 'LPP-WMMSE', 'WMMSE');
xlabel('SNR (dB)');
ylabel('系统吞吐率（比特每用户）');
grid on;

sc    = [10^-1, 10^-2, 10^-3, 10^-4, 10^-5, 10^-6, 10^-7];
sbcd  = [1, 8, 25, 59, 125, 203, 307];
pbcd6 = [1, 4, 12, 28,  48,  79, 119];
pbcd4 = [1, 3,  9, 22,  37,  56,  82];
semilogx(sc, sbcd, '-^', sc, pbcd6, '-o', sc, pbcd4, '-s', 'LineWidth', 1.5);
legend('连续BCD', '并行BCD, \beta_{ik}=6', '并行BCD, \beta_{ik}=4');
xlabel('停止条件');
ylabel('迭代次数');
grid on;

