function [norme, taux] = tauxL2(u_true, u_h, x, h_max, alpha)
	C = L2(u_true - u_h, x);

	norme= C * h_max^alpha;
	taux= -1 * log(C) - alpha*log(1/h_max);
end
