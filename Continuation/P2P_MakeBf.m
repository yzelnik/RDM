function [bifdata]=P2P_MakeBf(foldername,Ps,Es,AnlFunc,varargin)
% Read many P2P states and create a bif list
% Vs=P2P_MakeBf(foldername,Ps,Es,AnfFunc)

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:});

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

%bifdata=(dirlist);
for fileind=1:length(nums)
	% Read state, and then run AnlFunc on it
	[st,pp] = P2P_ReadSt(sprintf('%s/%s%d',foldername,basename,nums(fileind)),Ps,Es);
	%disp(sprintf('%s/%s%d%s',foldername,basename,nums(fileind),'.mat'));
	prebf1=bradat(pp);	
	prebf2=pp.fuha.outfu(pp,pp.u);
	Ps.mu=pp.u(end-4)+1;
	%disp(Ps.mu)
	[bf1,bf2] = AnlFunc(st,Ps,Es);
	bifdata(fileind,:)=[prebf1(:)' prebf2(:)' bf2(:)' bf1];
	%Ps
end;



end
