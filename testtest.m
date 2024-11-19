clear;clc
S = RandStream('mt19937ar','Seed',5489);
sigin = sqrt(2)*sin(0:pi/8:6*pi);
sigout1 = awgn(sigin,10,0,S);
sigout2 = awgn(sigin,10,0,S);
isequal(sigout1,sigout2)
reset(S);
sigout3 = awgn(sigin,10,0,S);
isequal(sigout1,sigout3)
plot(sigin)