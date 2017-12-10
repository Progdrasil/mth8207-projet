function [norme] = L2(u, x)
	norme = sqrt(trapz(x, u .^ 2));
end
