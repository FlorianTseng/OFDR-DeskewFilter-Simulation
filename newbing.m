
% 参数
N = 1000; % 采样点数
f0 = 1e9; % 载波频率
B = 10e6; % 带宽
T = N / B; % 扫描时间
t = linspace(0, T, N); % 时间向量
f = linspace(-B/2, B/2, N); % 频率向量

% 信号生成
s = exp(1j * (2 * pi * f0 * t + 20 * sin(2 * pi * 1e6 * t)));

% 去斜滤波器设计
H = exp(-1j * pi * B * t.^2 / T);

% 去斜滤波
s_out = ifft(fft(s) .* H);

% 绘图
figure;
subplot(2,2,1);
plot(t, real(s));
title('原始信号');
xlabel('时间 (s)');
ylabel('幅度');

subplot(2,2,2);
plot(f, abs(fftshift(fft(s))));
title('原始信号频谱');
xlabel('频率 (Hz)');
ylabel('幅度');

subplot(2,2,3);
plot(t, real(s_out));
title('去斜后信号');
xlabel('时间 (s)');
ylabel('幅度');

subplot(2,2,4);
plot(f, abs(fftshift(fft(s_out))));
title('去斜后信号频谱');
xlabel('频率 (Hz)');
ylabel('幅度');