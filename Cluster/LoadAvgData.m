function [data,valnum] = LoadAvgData(filename,column)
if(nargin<2) column=1; end;

load(filename);
if(isfield(Es,'PartsCollected'))
    TotBf = TotBf(logical(Es.PartsCollected),:);
end;
collen=size(TotBf,2);
if(length(column)==1)
    axdata=unique(TotBf(:,column));
    data = zeros(length(axdata),size(TotBf,2));
    valnum=zeros(size(axdata));
    for ii=1:length(axdata)
        tmp = TotBf(TotBf(:,column)==axdata(ii),:);
        valnum(ii) = size(tmp,1); % how many did we get?
        for kk=1:collen
        	tmp3 = tmp(:,kk);
        	tmp3 = tmp3(~isinf(tmp3)&~isnan(tmp3)); % make sure no nan's and inf's
        	data(ii,kk) = mean(tmp3); % calculate average
        end;
        
        %tmp = tmp(~isinf(tmp)&~isnan(tmp)); % make sure no nan's and inf's
        %data(ii,:)=mean(tmp,1); % calculate average
    end;
elseif(length(column)==2)
    ax1=unique(TotBf(:,column(1)));
    ax2=unique(TotBf(:,column(2)));
    data = zeros(length(ax1),length(ax2),size(TotBf,2));
    valnum=zeros(length(ax1),length(ax2));
        for ii=1:length(ax1) % go over axis 1
            tmp1 = TotBf(TotBf(:,column(1))==ax1(ii),:);
            for jj=1:length(ax2) % go over axis 2
                tmp2=tmp1(tmp1(:,column(2))==ax2(jj),:);
                valnum(ii,jj) = size(tmp2,1); % how many did we get?
                for kk=1:collen
                    tmp3 = tmp2(:,kk);
                    tmp3 = tmp3(~isinf(tmp3)&~isnan(tmp3)); % make sure no nan's and inf's
                    data(ii,jj,kk) = mean(tmp3); % calculate average
                end;
            end;
        end;
else
    error('only 1 or 2 columns supported');
end;

end
