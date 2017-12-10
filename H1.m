function [norme] = H1(u, x, h)
	up = diff(u)/h;
	nUp = L2(up, x(1:length(up)));
	nU = L2(u, x);
	norme = sqrt(nU^2 + nUp^2);
end
