function [newVect] = concat(varargin)
	k = 0;
	for i = 1:length(varargin)
		tempVect = varargin{i};
		for j = 1:length(tempVect)
			k = k+1;
			newVect(k) = tempVect(j);
		end
	end
end
