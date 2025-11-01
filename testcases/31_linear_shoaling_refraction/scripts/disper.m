function [k,C,Cg]=disper(h,T)
% [k,C,Cg]=disper(h,T)
% h - water depth
% T - wave period
% k - wave number
% C - wave celerity
% Cg - group velocity
% Approximate solution according to Guo (2002)
	g=9.81;
	sigma=2*pi./T;
	k = sigma.^2/g*(1-exp(-(sigma.*sqrt(h/g)).^2.5)).^(-0.4);
	C= sigma./k;
	n = 0.5+k.*h./sinh(2.*k.*h);
        Cg=n.*C;
