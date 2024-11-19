close all;clear all;clc

Z=10;n=0.9;c=3e5;%Z的单位为“千米”，千万不能用“米”！
tau=2*Z*n/c

f0 = 2.7619e6;              % 中心频率 f0 = 20/N*Fs=2.7619e6;
T = tau;                % 信号持续时间
B = 5.8e6;                  % 信号带宽
gamma = B/T;                    % 调频率
ratio = 5;                  % 过采样率
Fs = ratio*B;               % 采样频率
dt = 1/Fs;                  % 采样间隔
N = ceil(T/dt);             % 采样点数
t = ((0:N-1)-N/2)/N*T;      % 时间轴
f = ((0:N-1)-N/2)/N*Fs;     % 频率轴

Z=10;n=0.9;c=3e5;%Z的单位为“千米”，千万不能用“米”！
tau=2*Z*n/c;fb=gamma*tau;E0=5;
%f0是激光源初始频率；gamma是可调谐激光器调谐速度，单位Hz/s；
%fb是拍频大小，E0是光波振幅
dt=t+tau;
Er=E0*exp(j*(2*pi*f0*dt+pi*gamma*dt.^2+phi(dt)));
dEr=E0*exp(j*(2*pi*f0*t+pi*gamma*t.^2+phi(t)));
%phi是在t时刻光源随机波动的光相位，Er是本振参考光的光场
r=0.8;alpha=0.3;
R=r*exp(-alpha*tau*c/n);
Es=sqrt(R)*dEr;
It=abs(Er+Es).^2;
It2=E0^2*(1+R+2*sqrt(R)*cos(2*pi*(f0*tau+fb*t+0.5*gamma*tau^2+phi(t)-phi(dt))));
It-It2

figure(1)
subplot(331);
plot(t*1e6,angle(Er)),xlabel('时间 - [s]'),ylabel('相位'),title('本振参考光光场相位');
subplot(332);
plot(t*1e6,real(Er)),xlabel('时间 - [s]'),ylabel('幅值'),title('本振参考光光场实部振幅');
subplot(333);
plot(t*1e6,imag(Er)),xlabel('时间 - [s]'),ylabel('幅值'),title('本振参考光光场虚部振幅');

subplot(334);
plot(t*1e6,angle(Es)),xlabel('时间 - [s]'),ylabel('相位'),title('测试光光场相位');
subplot(335);
plot(t*1e6,real(Es)),xlabel('时间 - [s]'),ylabel('幅值'),title('测试光光场实部振幅');
subplot(336);
plot(t*1e6,imag(Es)),xlabel('时间 - [s]'),ylabel('幅值'),title('测试光光场虚部振幅');

subplot(337);
plot(t*1e6,angle(It)),xlabel('时间 - [s]'),ylabel('相位'),title('拍频信号相位');
subplot(338);
plot(t*1e6,real(It)),xlabel('时间 - [s]'),ylabel('幅值'),title('拍频信号实部振幅');
subplot(339);
plot(t*1e6,imag(It)),xlabel('时间 - [s]'),ylabel('幅值'),title('拍频信号虚部振幅');
