function []=Er(t,phi)
Er=E0*cos(2*pi*f0*t+pi*gamma*t.^2+phi);
end