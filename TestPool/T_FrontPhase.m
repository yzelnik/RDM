function [phase,share]=T_FrontPhase(Vs,Ps,Es,varargin)
% Check the relative phase of 2 variables in a front
% phase=T_FrontPhase(Vs,Ps,Es)
    
% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if( (~isfield(Es,'Vind')) || (Es.Vind<2) )
	Es.Vind = [1 2];
end;

[~,minvals,maxvals] = T_MinMax(Vs);                    % Get min & max values
middle1 = (minvals(Es.Vind(1))+maxvals(Es.Vind(1)))/2; % Value at mid-front for 1st variable
[~,midloc] = min(abs(Vs(:,Es.Vind(1))-middle1));       % Where is it?

tmp=abs(Vs(:,Es.Vind(1))-middle1);
[val1,midloc1] = min(tmp);
tmp(midloc1)=inf;
[val2,midloc2] = min(tmp);
share1 = val2/(val1+val2);
share2 = val1/(val1+val2);
valatmid = Vs(midloc1,Es.Vind(2))*share1+Vs(midloc2,Es.Vind(2))*share2;
%[midloc1 midloc2 share1 share2 val1*share1+val2*share2 valatmid]
%valatmid = Vs(midloc,Es.Vind(2));
share = [share1 share2];

middle2 = (minvals(Es.Vind(2))+maxvals(Es.Vind(2)))/2; % Value at mid-front for 2nd variable
phase(1) = 2*(valatmid-middle2)/(maxvals(Es.Vind(2))-minvals(Es.Vind(2)));
%[middle1 middle2 midloc Vs(midloc,Es.Vind(2)) [(Vs(midloc,Es.Vind(2))-middle2) (maxvals(Es.Vind(2))-minvals(Es.Vind(2)))]]

% And now just considering the ends
minvals=Vs(1,:,1);
maxvals=Vs(end,:,1);
middle1 = (minvals(Es.Vind(1))+maxvals(Es.Vind(1)))/2; % Value at mid-front for 1st variable
[~,midloc] = min(abs(Vs(:,Es.Vind(1))-middle1));       % Where is it?

tmp=abs(Vs(:,Es.Vind(1))-middle1);
[val1,midloc1] = min(tmp);
tmp(midloc1)=inf;
[val2,midloc2] = min(tmp);
share1 = val2/(val1+val2);
share2 = val1/(val1+val2);
%share1=max(min(share1,0.7),0.3);
%share2=1-share1;

valatmid = Vs(midloc1,Es.Vind(2))*share1+Vs(midloc2,Es.Vind(2))*share2;
%disp([midloc1 midloc2 share1 share2 val1 val2 val1*share1+val2*share2 Vs(midloc1,Es.Vind(2)) Vs(midloc2,Es.Vind(2)) valatmid])
%valatmid = Vs(midloc,Es.Vind(2));

middle2 = (minvals(Es.Vind(2))+maxvals(Es.Vind(2)))/2; % Value at mid-front for 2nd variable
phase(2) = 2*(valatmid-middle2)/(maxvals(Es.Vind(2))-minvals(Es.Vind(2)));
%disp([middle1 middle2 midloc Vs(midloc,Es.Vind(2)) [(Vs(midloc,Es.Vind(2))-middle2) (maxvals(Es.Vind(2))-minvals(Es.Vind(2)))]])

share = [share share1 share2];

end

