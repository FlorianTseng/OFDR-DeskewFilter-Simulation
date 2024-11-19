function [y] = stretch(nscat,taup,f0,b,scat_range,rrec,scat_rcs,winid)
eps = 1e-16;
htau = taup/2;
c = 3e8;
trec = 2.*rrec/c;
n = fix(2.*trec*b);%fix向0取整
m = power_integer_2(n);
nfft = 2.^m;
%初始化矩阵
x(nscat,1:n) = 0;
y(1:n) = 0;
%确定适当的窗口
if (winid == 0)
    win(1:n) = 1;
    win = win';    %不加窗
else
    if(winid == 1)
    win = hamming(n);   %加汉明
    else
        if(winid == 2)
            win = kaiser(n,pi);  %参数为pi的kaiser窗
        else
            if(winid == 3)
                win = chebwin(n,60);   %旁瓣为-60db的切比雪夫窗
            end
        end
    end
end

deltar = c/2/b;
max_rrec = deltar*nfft/2;
maxr = max(scat_range);

if(rrec>max_rrec|maxr>=rrec)
    'Error.Receive window is too large;or scatterers fall outside window';
    return
end

t = linspace(0,taup,n);
for j = 1:1:nscat
    range = scat_range(j);
    psi1 = 4.*pi*range*f0/c-4.*pi*b*range*range/c/c/taup;
    psi2 = (2*4.*pi*b*range/c/taup).*t;
    x(j,:) = scat_rcs(j).*exp(i*psi1+i.*psi2);
    y = y + x(j,:);
end

figure
plot(t,real(y),'k');
xlabel('时间延迟');
ylabel('未压缩的回波信号');

ywin = y.*win';
yfft = fft(y,n)./n;
out = fftshift(abs(yfft));
figure
delinc = rrec/n;
dist = linspace((-rrec/2),rrec/2,n);
plot(dist,out,'k');
xlabel('相对距离/m');
ylabel('压缩回波信号');

