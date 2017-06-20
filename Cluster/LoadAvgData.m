function [data,valnum] = LoadAvgData(filename,column)
if(nargin<2) column=1; end;

load(filename);
if(isfield(Es,'PartsCollected'))
    TotBf = TotBf(logical(Es.PartsCollected),:);
end;

if(length(column)==1)
    axdata=unique(TotBf(:,column));
    data = zeros(length(axdata),size(TotBf,2));
    valnum=data;
    for ii=1:length(axdata)
        tmp = TotBf(TotBf(:,column)==axdata(ii),:);
        valnum(ii,1) =size(tmp,1);
        data(ii,:)=mean(tmp,1);
    end;
elseif(length(column)==2)
    ax1=unique(TotBf(:,column(1)));
    ax2=unique(TotBf(:,column(2)));
    data = zeros(length(ax1),length(ax2),size(TotBf,2));
    valnum=data;
        for ii=1:length(ax1) % go over axis 1
            tmp1 = TotBf(TotBf(:,column(1))==ax1(ii),:);
            for jj=1:length(ax2) % go over axis 2
                tmp2=tmp1(tmp1(:,column(2))==ax2(jj),:);
                data(ii,jj,:)=mean(tmp2,1); % calculate average
                valnum(ii,jj)=size(tmp2,1);
            end;
        end;
else
    error('only 1 or 2 columns supported');
end;

end
