function taux = tauxConv(e_pres, e_last, h_pres, h_last)
%	tauxConv function
%	inputs:
%		- e_pres = norme of error of present analysis
%		- e_last = norme of error of last analysis
%		- h_pres = max element size of present analysis
%		- h_last = max element size of last analysis
	taux= log(e_last / e_pres) / log(h_last / h_pres);
end
