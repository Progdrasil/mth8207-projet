function [norme, taux] = tauxConv(u_true, u_h, x, h_max, alpha, type_norme)
	type_norme = upper(type_norme);
	if (~strcmp(type_norme, 'L2') && ~strcmp(type_norme, 'H1'))
		error('Type de norme pas reconnaissable');
	end
	fh = str2func(type_norme);
	if (strcmp(type_norme, 'H1'))
		C = fh(u_true - u_h, x, h_max);
	else
		C = fh(u_true - u_h, x);
	end

	norme= C * h_max^alpha;
	taux= -1 * log(C) - alpha*log(1/h_max);
end
