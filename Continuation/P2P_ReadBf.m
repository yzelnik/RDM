function [bifdata]=P2P_ReadBf(foldername,Ps,Es,varargin)
% Read P2P files and read their bif data
% bifdata=P2P_ReadBf(filename,Ps,Es)

% Update online if necessary
[~,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:});


bifdata=[];
basename = 'p';

%ffname=[pre '/' fname '.mat'];
fullpath = [foldername '/' basename '*'];
dirlist=dir(fullpath);
if(isempty(dirlist))
	dirlist=dir(foldername);
	if(isempty(dirlist))
		error(sprintf('Could not find state files using path "%s"',fullpath));
	else
		warning(sprintf('Did not find files in full path, instead used: "%s"',foldername));
	end;
end;
nums=[];
for fileind=1:length(dirlist)

	[aa,bb]=strread(dirlist(fileind).name,'%[^0123456789]%d.mat');
	if(~isempty(bb))
		nums=[nums bb];
	end;
end;
nums=sort(nums);

[st,pp] = P2P_ReadSt(sprintf('%s/%s%d',foldername,basename,nums(end)),Ps,Es);
bifdata=pp.branch';




end





%bifdata=[];
%basename = 'p';


%if(strcmp(filename(end-3:end),'.mat')
%    [st,pp] = P2P_ReadSt(filename(1:end-4),Ps,Es);
%else
%dirlist=dir(filename);

%if(isempty(dirlist))
%    foldername=filename;
%disp(1)
%end;
%ffname=[pre '/' fname '.mat'];
%fullpath = [foldername '/' basename '*'];
%dirlist=dir(fullpath);
%if(isempty(dirlist))%
%	dirlist=dir(foldername);
%	if(isempty(dirlist))
%		error(sprintf('Could not find state files using path "%s"',fullpath));
%	else
%		warning(sprintf('Did not find files in full path, instead used: "%s"',foldername));
%	end;
%end;
%nums=[];
%for fileind=1:length(dirlist)
%
%	[aa,bb]=strread(dirlist(fileind).name,'%[^0123456789]%d.mat');
%	if(~isempty(bb))
%		nums=[nums bb];
%	end;
%end;
%nums=sort(nums);

