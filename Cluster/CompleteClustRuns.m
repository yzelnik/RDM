function CompleteClustRuns(partruns,Vs)
% get partial runs data, from a structure or file
if(isstruct(partruns))
    tmpdata=partruns;
elseif(isstring(partruns) || ischar(partruns))
    tmpdata=load(partruns);
else
    error('Either a filename or strcture is needed.');
end;

if(~isfield(tmpdata,'Vs')) % If no Vs exists in tmpdata
    if(nargin>1)
        tmpdata.Vs=Vs;
    else
        tmpdata.Vs=0;
        warning('No Vs variables was given, using Vs=0 as default');
    end;
end;

if(~isfield(tmpdata.Es,'PartsCollected'))
    error('Es.PartsCollected is not defined. Did all runs end successfully?');
end

if(length(tmpdata.Es.PartsCollected)<size(tmpdata.TotBf,1))
    tmpdata.Es.PartsCollected(size(tmpdata.TotBf,1))=0;
end;

% indices of missing parts
missingparts = find(1-tmpdata.Es.PartsCollected);

disp(sprintf('found %d missing parts out of %d, starting to run.',length(missingparts),size(tmpdata.TotBf,1)));

tmpdata.Es.RunsChoice = missingparts;

if(isempty(tmpdata.TotSt))
    [~,missingbf] = runpar(tmpdata.Vs,tmpdata.Ps,tmpdata.Es,'Es.Verbose',1);
else
    [missingst,missingbf] = runpar(tmpdata.Vs,tmpdata.Ps,tmpdata.Es,'Es.Verbose',1);
    tmpdata.TotSt(missingparts)=missingst;
end;

tmpdata.TotBf(missingparts,:)=missingbf;

if(isstring(partruns) || ischar(partruns))
    basename = partruns;
    if(length(basename)>4) && strcmp(basename(end-3:end),'.mat')
        basename = basename(1:end-4);
    end;
else
    basename = sprintf('ClustData%04d',randi(1e4)-1);
end;

Ps=tmpdata.Ps;
Es=tmpdata.Es;
TotBf=tmpdata.TotBf;
TotSt=tmpdata.TotSt;

save(sprintf('%s_completed',basename),'Ps','Es','TotBf','TotSt');

end
