function [varargout]=ExtractMultFile(basename,ranges,vars)
% use: [var1,var2,...]=ExtractMultFile(basename,ranges,vars)
% basename should look like: run_%d_prm%d, with %d for each range
% vars contains the names of variables to read into var1,... data sets

% Force into cell structure
if(~iscell(ranges)) ranges={ranges}; end;
if(~iscell(vars)) vars={vars}; end;

% Figure out the total number of files to load
totsize=1; 
for ii=1:length(ranges) 
    basesize(ii)=length(ranges{ii}); 
    totsize=totsize*basesize(ii); 
end;
if(length(basesize)==1)
    basesize = [basesize 1];
end;
% Make an organized list of all file indices
tmpsize=ones(size(basesize));
for ii=1:length(ranges) 
    basesize(ii)=1; 
    tmpsize(ii)=length(ranges{ii}); 
    tmp=repmat(reshape(ranges{ii},tmpsize),basesize); 
    finind(:,ii)=reshape(tmp,totsize,1); 
    basesize(ii)=length(ranges{ii}); 
    tmpsize(ii)=1; 
end;

count=0;
% Go over list and read files
for ii=1:size(finind,1)
    %disp(sprintf(basename,finind(ii,:)));
    fullname = sprintf(basename,finind(ii,:));
    if(exist(fullname,'file'))
        count=count+1;
    	tmpload=load(sprintf(basename,finind(ii,:)));
        for jj=1:length(vars)
            varargout{jj}{ii}=tmpload.(vars{jj});
        end;
    else
        for jj=1:length(vars)
            varargout{jj}{ii}=[];
        end;
    end;
end;
disp(sprintf('Found %d out of %d files.',count,size(finind,1)));

% Reshape final result
for jj=1:length(vars)
        varargout{jj}=reshape(varargout{jj},basesize);
end;

end
