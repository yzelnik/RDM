function [Out,Score]=CheckSS(Vs,Ps,Es,varargin)
% Check convergence to steady state
% Out=0: no convergence
% Out=1: convergence

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

V1 = Vs(:,:,1);
V2 = Vs(:,:,2);

Score = max(max(abs((V2-V1))));%./max(abs(V2),Es.StSmall))));
%disp([max(abs(V2(1)),Es.StSmall)])% Score mean(mean(abs((V2-V1)./max(abs(V2),Es.StSmall))))])
if Score<Es.SsThresh
    Out=1;
else
    Out=0;
end

end

