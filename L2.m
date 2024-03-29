function [norme] = L2(uh, xh, ut, xt, nEls, xMax, xMin)
	if (nargin == 2)
		norme = sqrt(trapz(xh, uh .^ 2));
		return;
	end

	% if (nargin < 5)
		nG = 3;
	% end

	% premier points
	k=1;
	xe(k) = xh(1);
	uhj = uh(1);
	utj = interp1(xt, ut, xe(k));
	e(k) = (utj - uhj) ^ 2;

	dxG = (xMax - xMin) / (nEls * nG);

	for i = 1: length(nEls)
		X(i) = xe(k);
		for j = 1:nG
			k = k + 1;
			xe(k) = X(i) + j * dxG;
			uhj = interp1(xh, uh, xe(k));
			utj = interp1(xt, ut, xe(k));
			e(k) = (utj - uhj) ^ 2;
			if isnan(e(k)) || isnan(xe(k))
				k = k-1;
				error('fuck this shit')
			end
		end
	end

	norme = sqrt(trapz(xe, e));
	return;
end
