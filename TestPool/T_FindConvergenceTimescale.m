function [tscale,tottime,totdiff] = T_FindConvergenceTimescale(Vs,Ps,Es,varargin)
% Run Integration of a numerical solution until convergence
% estimate the time scale (assuming exponential decay)
% Also return the total time of integration, and the difference between initial and final states

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

% Initialization
Es.Tmax  = Es.Tdest;

[~,his,tottime]=run2ss(Vs,Ps,Es,'Es.OLdraw',1,'Es.TestFunc',@T_L2Norm);

totdiff  = diff(his([1 end],3));

[~,mind] = max(his(:,3));		% Find the max difference over time, assume convergence starts after it	
coeffs   = polyfit(his(mind:end,1),log(his(mind:end,3)),1);
tscale   = -1/coeffs(1);		% use the slope (of a logarithmic function to begin with) to get the timescale
plot(his(mind:end,1),log(his(mind:end,3)))

end