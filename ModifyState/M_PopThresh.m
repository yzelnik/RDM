function VsOut = M_PopThresh(Vs,Ps,Es,varargin)
% Cut a given variable down to zero if it falls below a threshold in the whole system
% VsOut = M_PopThresh(Vs,Ps,Es)
% Threshold is set by Es.PopThresh, (default is set by Es.StSmall)
% If Es.PopThresh is a vector, it is done per variable, otherwise for Es.VarInd

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Es=InsertDefaultValues(Es,'PopThresh',Es.StSmall);

VsOut=Vs;

if(length(Es.PopThresh)==1)
    %if(max(VsOut(:,Es.VarInd,1))<Es.PopThresh)
    %    VsOut(:,Es.VarInd,1)=0;
    %end;
    sites=VsOut(:,Es.VarInd,1)<Es.PopThresh;
    VsOut(sites,Es.VarInd,1)=0;
else
    for ii=1:length(Es.PopThresh)
    %    if(max(VsOut(:,ii,1))<Es.PopThresh(ii))
    %        VsOut(:,ii,1)=0;
    %    end;
        sites=VsOut(:,ii,1)<Es.PopThresh(ii);
        VsOut(sites,ii,1)=0;
    end;
end;


end
