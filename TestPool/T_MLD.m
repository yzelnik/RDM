function [res,ratio]=T_MLD(Vs,Ps,Es,varargin)
% Check which (2nd or 3rd) state is More Linear Dependent on the 1st state 
% [res,ratio]=T_MLD(Vs,Ps,Es)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;


if( (~isfield(Es,'VarInd')) | (Es.VarInd<3) )
	Es.VarInd = [1 2 3];
end;

dif12 = sum(abs(normcol(Vs(:,Es.VarInd(1),1))-normcol(Vs(:,Es.VarInd(2),1))));
dif13 = sum(abs(normcol(Vs(:,Es.VarInd(1),1))-normcol(Vs(:,Es.VarInd(3),1))));

ratio = log10(dif12/dif13);
res = ratio>0;

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


