clc;close all;
c=3e8;n=0.9;dv=(0:100)/100;
z=c/pi/n./dv/1e6;
plot(dv,z);
title('光源线宽与OFDR系统最长测量距离的关系'),xlabel('光源线宽 - [MHz]'),ylabel('系统最长测量距离 - [m]');