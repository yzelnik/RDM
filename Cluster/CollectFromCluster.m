function CollectFromCluster(basefilename,subfilenum,varargin)
% load a set of files and merge them into one
if(nargin<2)
    subfilenum=0;
end;
PartsCollected=[]; % keep track of things

FileType='*';
if(length(basefilename)>4) % If a predefined file type is used 
    if(strcmp(basefilename(end-3:end),'.mat') || strcmp(basefilename(end-3:end),'.csv') || strcmp(basefilename(end-3:end),'.txt'))
        FileType=basefilename(end-3:end);
        basefilename=basefilename(1:end-4);    
    end;
end;
if(subfilenum)
    NumTxt=num2str(subfilenum);
else
    NumTxt='*';
end;

list=dir(sprintf('*%s_part*of%s%s',basefilename,NumTxt,FileType));

TotSt={};
TotBf=[];
if(isempty(list))
    error('Could not find files matching: *%s_part*of%s%s',basefilename,NumTxt,FileType);
else
    if(strcmp(NumTxt,'*')||strcmp(FileType,'*'))
        filename = list(1).name;
        FileType=filename(end-3:end);
        tmpnuminds=regexp(filename,'[0-9]');
        tmpend=find(diff(tmpnuminds)~=1);
        NumTxt=filename(tmpnuminds(tmpend(end)+1:end));
        list=dir(sprintf('*%s_part*of%s%s',basefilename,NumTxt,FileType));
    end;
   
    if(strcmp(FileType,'.mat'))
        for ii=1:length(list)
            tmp=load(list(ii).name);
        
            %inds = tmp.Es.RunsChoice;
            Es=SortOutBfParameters(tmp.Es);
            totvals = Es.BfVal; % list of parameter-combinations
            inds = ceil(size(totvals,1)*(Es.RunsChoice(1)-1)/Es.RunsChoice(2))+1 : ceil(size(totvals,1)*Es.RunsChoice(1)/Es.RunsChoice(2));
            if(isfield(tmp,'StData')) 
                datalen=min(size(tmp.StData,1),size(tmp.BfData,1));
            else
                datalen=size(tmp.BfData,1);
            end;
            
            if(datalen<length(inds))
                warning('Could not load all data from %s, only %d out of %d run-results available.',list(ii).name,datalen,length(inds));
                %if(isempty(PartsCollected))
                %    PartsCollected(1:size(TotBf,1))=1; % first "access"
                %end;
                %PartsCollected(size(TotBf,1)+(1:datalen))=1;
            end;
            PartsCollected(inds(1:datalen))=1; % mark parts are correctly read
            if(isfield(tmp,'StData') && (~isempty(tmp.StData)))
                TotSt(inds(1:datalen))=tmp.StData(1:datalen);
            end;
            TotBf(inds(1:datalen),1:size(tmp.BfData,2))=tmp.BfData(1:datalen,:);
        end;
        Ps = tmp.Ps;
        Es = tmp.Es;
        if(sum(PartsCollected)<size(TotBf,1)) 
            Es.PartsCollected=PartsCollected;   % Used if there were missing parts
        end;
        save(sprintf('%s_totalof%d%s',basefilename,length(list),FileType),'TotSt','TotBf','Es','Ps');
    else
        totbf=[];
        for ii=1:length(list)
            data = dlmread(list(ii).name);
            totbf(size(totbf,1)+(1:size(data,1)),1:size(data,2)) = data; 
        end;
        dlmwrite(sprintf('%s_totalof%d.%s',basefilename,length(list),FileType),totbf);
    end;
    
end;

end

