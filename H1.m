function [norme] = H1(uh, Duh, xh, ut, Dut, xt, nEls, xMax, xMin)

	% Duh = diff(uh)/h;

	if (nargin <= 3)
		% calculs normes L2
		nUp = L2(Duh, xh(1:length(Duh)));
		nU = L2(uh, xh);
		% calculs norme H1
		norme = sqrt(nU^2 + nUp^2);
		return;
	else
		% if (nargin < 7)
			nG = 3;
		% end

		% Dut = diff(ut)/h;
		[Dut, Duh] = normalise(Dut, Duh, false);

		% calculs normes L2
		nEp = L2(Duh, xh(1:length(Duh)), Dut, xt(1:length(Dut)),  nEls, xMax, xMin);
		nE = L2(uh, xh, ut, xt,  nEls, xMax, xMin);

		% calcul norme H1
		norme = sqrt(nE^2 + nEp^2);

		return;
	end
end
