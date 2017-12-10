function [norme] = H1(u, x, h)
	up = gradient(u);
	disp(up)
	disp(L2(up, x(1:length(up))))
	norme = sqrt((L2(u, x) ^ 2) + (L2(up, x(1:length(up))) ^ 2));
end
