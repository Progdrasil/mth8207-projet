function [norme] = L2(uh, xh, ut, xt, nG)
	if (nargin == 2)
		norme = sqrt(trapz(xh, uh .^ 2));
		return;
	end

	if (nargin < 5)
		nG = 3;
	end

	% premier points
	k=1;
	xe(k) = xh(1);
	uhj = uh(1);
	utj = interp1(xt, ut, xe(k));
	e(k) = (utj - uhj) ^ 2;

	for i = 2: length(xh)
		dxG = (xh(i) - xh(i-1)) / 3;
		for j = 1:nG
			k = k + 1;
			xe(k) = xh(i-1) + j * dxG;
			uhj = interp1(xh, uh, xe(k));
			utj = interp1(xt, ut, xe(k));
			e(k) = (utj - uhj) ^ 2;
			if isnan(e(k))
				k = k-1;
			end
		end
	end

	norme = sqrt(trapz(xe, e));
	return;
end
