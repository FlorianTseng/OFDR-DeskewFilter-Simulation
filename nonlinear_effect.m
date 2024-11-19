Itideal=E0^2*(1+R+2*sqrt(R)*cos(2*pi*(f0*tau+fb*t+0.5*gamma*tau^2)));
spec=fft(It);
spec=fftshift(spec);
plot(t/n*c,abs(spec),'b')
hold on
spec=fft(Itideal);
spec=fftshift(spec);
plot(t/n*c,abs(spec),'r');
xlabel('距离 - [m]');ylabel('振幅');
legend('非线性调谐效应','理想条件');