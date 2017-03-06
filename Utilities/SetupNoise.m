function [noisefunc,noiseparm]=SetupNoise(Ps,Es)
% A Utility function to setup noise for integration
% [noisefunc,noiseparm]=SetupNoise(Ps,Es)

if(isempty(Ps.NoiseType))
    error('No noise type defined by Ps.NoiseType');
end;

% decide on noise function
if(isfloat(Ps.NoiseType))
    % simple power noise function (env noise y(2)=1, demographic noise y(2)=0.5, etc)
    noisefunc = @(x,y) y(1)*x.^y(2); 
    noiseparm = Ps.NoiseType;
elseif(iscell(Ps.NoiseType))
    noisefunc = Ps.NoiseType{1};
    noiseparm = Ps.NoiseType{2};
else % assume it is a function
    noisefunc = Ps.NoiseType;
    noiseparm = 0;
end;   

% make sure parameters are set up correctly for noise
if(~(size(noiseparm,2)==Ps.VarNum))
	if(size(noiseparm,1)==Ps.VarNum)
        noiseparm=noiseparm'; % flip column to row
    else
        noiseparm=repmat(noiseparm(:,1),1,Ps.VarNum);
	end;
end;

% pad parameters with zeros
noiseparm = [noiseparm ; zeros(1,size(noiseparm,2))];


end
