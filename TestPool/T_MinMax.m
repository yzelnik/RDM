function [minmaxvals,minvals,maxvals,avgvals] = T_MinMax(Vs,varargin)
% Calculate the Min and Max values of the state variables
% [minmaxvals,minvals,maxvals,avgvals] = T_MinMax(Vs)

for ii=1:size(Vs,2)
	minvals(ii) = min(Vs(:,ii,1));
	maxvals(ii) = max(Vs(:,ii,1));
	minmaxvals(ii) = maxvals(ii)-minvals(ii);
	avgvals(ii) = mean(Vs(:,ii,1));
end;

end