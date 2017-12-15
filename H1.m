function [norme] = H1(uh, Duh, xh, h, ut, xt, nG)

	% Duh = diff(uh)/h;

	if (nargin <= 4)
		% calculs normes L2
		nUp = L2(Duh, xh(1:length(Duh)));
		nU = L2(uh, xh);
		% calculs norme H1
		norme = sqrt(nU^2 + nUp^2);
		return;
	else
		if (nargin < 7)
			nG = 3;
		end

		Dut = diff(ut)/h;

		% calculs normes L2
		nEp = L2(Duh, xh(1:length(Duh)), Dut, xt(1:length(Dut)), nG);
		nE = L2(uh, xh, ut, xt, nG);

		% calcul norme H1
		norme = sqrt(nE^2 + nEp^2);

		return;
	end
end
