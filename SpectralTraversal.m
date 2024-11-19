clear all;close all;clc
zmin=0;
length=1000;
co1=zeros(1,length);
co2=co1;
for Z=zmin:(zmin+length)
n=0.9;c=3e8;     %Z为光纤测量目标点距离
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
zf=0.0001;pl=B;

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

%利用希尔伯特变换，将拍频信号I(t)转换到复指数形式
It=2*sqrt(R)*E0^2*exp(1i*2*pi*(f0*tau+fb*t+0.5*gamma*tau.^2+phi(t,zf,pl)-phi(dt,zf,pl)));
%It1=It.*conj(se(t));   
It1=2*sqrt(R)*E0^2*exp(1i*2*pi*(f0*tau+fb*t+0.5*gamma*tau.^2-phi(dt,zf,pl)));

It2=ifft(fft(It1).*exp(1i*pi*f.^2/gamma));
%It2=2*sqrt(R)*E0^2*exp(j*2*pi*(f0*tau+fb*t))*s(t);
s=ifft(fft(conj(se(t,zf,pl))).*exp(1i*pi*f.^2/gamma));                     %信号卷积去斜滤波器
It3=It2.*conj(s);

subplot(121)

spectrum=fft(It3);
spec=spectrum.*exp(1i*pi*f.^2/gamma);
spec=fftshift(spec);
plot(t/n*c,abs(spec),'r');
xlabel('距离 - [m]');ylabel('振幅');
hold on
spec=fft(It);
spec=fftshift(spec);
plot(t/n*c,abs(spec),'b')
legend('去斜滤波器','原始频谱');
title(sprintf('相位噪声振幅为%d时拍频信号散射峰谱',zf));
axis([zmin zmin+length 0 inf]);

coematrix=corrcoef(t,phase(It));
co1(Z-zmin+1)=abs(coematrix(1,2));
coematrix=corrcoef(t,phase(It3));
co2(Z-zmin+1)=abs(coematrix(1,2));

end
subplot(122)
x=zmin:(zmin+length);
%p1=polyfit(x,co1(1,x-zmin+1),9);
plot(x,co1(1,x-zmin+1),'*b');
hold on
plot(x,co2(1,x-zmin+1),'+r');
hold on
plot(x,co1(1,x-zmin+1));
%y1=polyval(p1,x);
hold on
%plot(x,y1);
hold on

%p2=polyfit(x,co2(1,x-zmin+1),9);
%plot(x,co2(1,x-zmin+1),'+r');
%y2=polyval(p2,x);
plot(x,co2(1,x-zmin+1));
hold on
%plot(x,y2);
xlabel('距离 - [m]');ylabel('相关系数');
axis([zmin zmin+length 0 inf]);
title('去斜滤波前后拍频信号相位相关系数');
legend('原始信号','去斜滤波器');
