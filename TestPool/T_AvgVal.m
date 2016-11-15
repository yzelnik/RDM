function avgvals = T_AvgVal(Vs,varargin)
% Calculate the average values of the state variables
% avgvals = T_AvgVal(Vs,varargin)

avgvals = squeeze(mean(Vs(:,:,1),1));

end