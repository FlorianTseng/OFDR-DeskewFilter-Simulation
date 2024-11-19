clear all;clc;close all;

Z=1901;n=0.9;c=3e8;     %Z为光纤测量目标点距离
tau=2*Z*n/c;

T=tau;                % 信号持续时间
B=1e6;                  % 扫频范围，很关键
gamma=B/T;                    % 调频率
ratio=5000;                  % 过采样率
Fs = ratio*B;               % 采样频率
dt = 1/Fs;                  % 采样间隔
N = ceil(T/dt);             % 采样点数
t = (0:N-1)/N*T;      % 时间轴
f = (0:N-1)/N*Fs;     % 频率轴
f0 = 20/N*Fs;
zf=0.001;pl=B;

fb=gamma*tau;E0=5;
%f0是激光源初始频率；gamma是可调谐激光器调谐速度，单位Hz/s；
%fb是拍频大小，E0是光波振幅
dt=t+tau;
Er=E0*exp(1i*(2*pi*f0*t+pi*gamma*t.^2+phi(t,zf,pl)));
dEr=E0*exp(1i*(2*pi*f0*dt+pi*gamma*dt.^2+phi(dt,zf,pl)));
%phi是在t时刻光源随机波动的光相位，Er是本振参考光的光场
r=0.8;alpha=0.3/1000;
R=r*exp(-alpha*tau*c/n);
Es=sqrt(R)*dEr;         %Es为测试光光场
It=E0^2*(1+R+2*sqrt(R)*cos(2*pi*(f0*tau+fb*t+0.5*gamma*tau^2+phi(t,zf,pl)-phi(dt,zf,pl))));

figure(1)
subplot(341);
plot(t*1e6,angle(Er)),xlabel('时间 - [\mus]'),ylabel('相位 - [rad]'),title('本振参考光光场相位');
subplot(342);
plot(t*1e6,real(Er)),xlabel('时间 - [\mus]'),ylabel('幅值'),title('本振参考光光场实部振幅');
subplot(343);
plot(t*1e6,imag(Er)),xlabel('时间 - [\mus]'),ylabel('幅值'),title('本振参考光光场虚部振幅');
subplot(344);
plot(t*1e6,2*pi*f0*dt+pi*gamma*dt.^2+phi(dt,zf,pl)),xlabel('时间 - [\mus]'),ylabel('相位 - [rad]'),title('本振参考光光场瞬时相位');
subplot(345);
plot(t*1e6,angle(Es)),xlabel('时间 - [\mus]'),ylabel('相位'),title('测试光光场相位');
subplot(346);
plot(t*1e6,real(Es)),xlabel('时间 - [\mus]'),ylabel('幅值'),title('测试光光场实部振幅');
subplot(347);
plot(t*1e6,imag(Es)),xlabel('时间 - [\mus]'),ylabel('幅值'),title('测试光光场虚部振幅');
subplot(348);
plot(t*1e6,2*pi*f0*t+pi*gamma*t.^2+phi(t,zf,pl)),xlabel('时间 - [\mus]'),ylabel('相位 - [rad]'),title('测试光光场瞬时相位');
subplot(349);
plot(t*1e6,angle(It)),xlabel('时间 - [\mus]'),ylabel('相位'),title('拍频信号相位');
subplot(3,4,10);
plot(t*1e6,real(It)),xlabel('时间 - [\mus]'),ylabel('幅值'),title('拍频信号实部振幅');
subplot(3,4,11);
plot(t*1e6,f0*tau+fb*t+0.5*gamma*tau^2+phi(t,zf,pl)-phi(dt,zf,pl)),xlabel('时间 - [\mus]'),ylabel('相位 - [rad]'),title('拍频信号瞬时相位');
subplot(3,4,12);
plot(t*1e6,2*pi*f0*dt+pi*gamma*dt.^2+phi(dt,zf,pl),'r');
hold on
plot(t*1e6,2*pi*f0*t+pi*gamma*t.^2+phi(t,zf,pl),'b');
xlabel('时间 - [\mus]'),ylabel('MHz'),title('LO和测试光频率曲线');
legend('LO','测试光');

%利用希尔伯特变换，将拍频信号I(t)转换到复指数形式
It=2*sqrt(R)*E0^2*exp(1i*2*pi*(f0*tau+fb*t+0.5*gamma*tau.^2+phi(t,zf,pl)-phi(dt,zf,pl)));
%It1=It.*conj(se(t));   
It1=2*sqrt(R)*E0^2*exp(1i*2*pi*(f0*tau+fb*t+0.5*gamma*tau.^2-phi(dt,zf,pl)));

It2=ifft(fft(It1).*exp(1i*pi*f.^2/gamma));
%It2=2*sqrt(R)*E0^2*exp(j*2*pi*(f0*tau+fb*t))*s(t);
s=ifft(fft(conj(se(t,zf,pl))).*exp(1i*pi*f.^2/gamma));                     %信号卷积去斜滤波器
It3=It2.*conj(s);
figure(2)
plot(t*1e6,phase(It3),'b');
hold on
plot(t*1e6,f0*tau+fb*t+0.5*gamma*tau^2+phi(t,zf,pl)-phi(dt,zf,pl),'r');
xlabel('时间 - [\mus]'),ylabel('相位 - [rad]'),title('拍频信号I(t)和I_{3}(t)的相位');
legend('I_{3}(t)相位','I(t)相位');

figure(3)
spectrum=fft(It3);
spec=spectrum.*exp(1i*pi*f.^2/gamma);
spec=fftshift(spec);
subplot(411);
plot(t/n*c,db(abs(spec)/max(spec)))
xlabel('距离 - [m]');ylabel('振幅 - [dB]');title(sprintf('目标距离为%dm时卷积去斜滤波器信号的瑞利散射峰距离谱（衰落）',Z));
grid on
subplot(412);
plot(t/n*c,abs(spec))
xlabel('距离 - [m]');ylabel('振幅');title(sprintf('目标距离为%dm时卷积去斜滤波器信号的瑞利散射散射峰距离谱（模值）',Z));

subplot(413)
spec=fft(It);
spec=fftshift(spec);
plot(t/n*c,db(abs(spec)/max(spec)));
xlabel('距离 - [m]');ylabel('振幅 - [dB]');title(sprintf('目标距离为%dm时拍频信号瑞利散射峰距离谱（衰落）',Z));
grid on
subplot(414)
plot(t/n*c,abs(spec))
xlabel('距离 - [m]');ylabel('振幅');title(sprintf('目标距离为%dm时拍频信号瑞利散射散射峰距离谱（模值）',Z));

figure(4)
spectrum=fft(It3);
spec=spectrum.*exp(1i*pi*f.^2/gamma);
spec=fftshift(spec);
plot(t/n*c,abs(spec),'r');
xlabel('距离 - [m]');ylabel('振幅');
hold on
spec=fft(It);
spec=fftshift(spec);
plot(t/n*c,abs(spec),'b')
legend('去斜滤波器','未经处理');
title(sprintf('距离%dm时拍频散射谱',Z));

corrcoef(t,phase(It))
corrcoef(t,phase(It3))
figure(5)
plot(t,phi(t,zf,pl)-phi(dt,zf,pl));
title(sprintf('距离%dm时相位噪声差值e(t)-e(t-τ)',Z));