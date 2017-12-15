function [newDphir, newDuEF] = normalise(dphir, duEF, fig)

	% derive fonction exacte
	normEX = max(abs(dphir));
	if abs(max(dphir)) < abs(min(dphir))
		normEX = -1 * normEX;
	end
	newDphir = dphir ./ normEX;

	% deriver fonction element fini
	normEF = max(abs(duEF));
	if abs(max(duEF)) < abs(min(duEF))
		normEF = -1 * normEF;
	end
	newDuEF = duEF ./ normEF;

	if fig
		figure
		hold on
		plot(r, newDphir)
		plot(xEF, newDuEF)
		hold off
		legend('EX', 'EF', 'location' , 'best')
	end

end
