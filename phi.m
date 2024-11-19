%非线性相位噪声
function [r]=phi(t,zf,pl)
r=zf*cos(2*pi*pl*t);
end