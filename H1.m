function [norme] = H1(uh, xh, h, ut, xt, nG)

	uhp = diff(uh)/h;

	if (nargin <= 3)
		% calculs normes L2
		nUp = L2(uhp, xh(1:length(uhp)));
		nU = L2(uh, xh);
		% calculs norme H1
		norme = sqrt(nU^2 + nUp^2);
		return;
	else
		if (nargin < 6)
			nG = 3;
		end

		utp = diff(ut)/h;

		% calculs normes L2
		nEp = L2(uhp, xh(1:length(uhp)), utp, xt(1:length(utp)), nG);
		nE = L2(uh, xh, ut, xt, nG);
		
		% calcul norme H1
		norme = sqrt(nE^2 + nEp^2);

		return;
	end
end
