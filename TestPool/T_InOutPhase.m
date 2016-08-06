function [res,sm]=T_InOutPhase(Vs,Ps,Es,varargin)
% Check if state is in or out of phase
% res=T_InOutPhase(Vs,Ps,Es)
% Returns 1 if the state is in phase, 0 otherwise

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;


if( (~isfield(Es,'Vind')) | (Es.Vind<2) )
	Es.Vind = [1 2];
end;

sm = mean((normcol(Vs(:,Es.Vind(1),1))-0.5).*(normcol(Vs(:,Es.Vind(2),1))-0.5)); 
res = sm>0;

end




function outmatt = normcol( matt )
%normalize values to [0 - 1] range, column by column
mn = min(matt);
mx = max(matt);
outmatt=zeros(size(matt));
for ii=1:size(matt,2);
    outmatt(:,ii)=(matt(:,ii)-mn(ii))/(mx(ii)-mn(ii));
end;

end


