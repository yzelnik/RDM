function [data,valnum] = LoadAvgData(filename,column)
if(nargin<2) column=1; end;

load(filename);
if(isfield(Es,'PartsCollected'))
    TotBf = TotBf(logical(Es.PartsCollected),:);
end;
axdata=unique(TotBf(:,column));

for ii=1:length(axdata)
    tmp = TotBf(TotBf(:,column)==axdata(ii),:);
    valnum(ii,1) =size(tmp,1);
    data(ii,:)=mean(tmp,1);
end;


end
