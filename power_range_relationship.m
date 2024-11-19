clc;close all;
r=50;
z=(0:r)/r;
a=0.3;k=0.5;x=0.99;h=a*1e-6;
c=sqrt(h*k*(1-k)*x);
s=c*exp(-a*z*r);
plot(z*r,s),xlabel('最长测试距离 - [km]'),ylabel('散射光功率 - [W]');
title('接收光功率与OFDR系统最长测量距离的关系');