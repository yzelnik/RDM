function [total,pervar] = T_L2Norm(Vs,varargin)
% Calculate the L2Norm of a given state

for ii=1:size(Vs,2)
	pervar(ii) = sqrt(mean(Vs(:,ii,1).^2));
end;
total = sqrt(sum(pervar.^2));

end