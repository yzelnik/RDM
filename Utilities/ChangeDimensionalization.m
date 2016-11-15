function [PsNew,VsNew,scaling]=ChangeDimensionalization(Vs,Ps,Es,NewParDef,varargin)
% Change a prameter setting to an alternative parameter setting
% For example, change from a dimentionalized to non-dimensionalized 
% [PsNew,VsNew,scaling]=ChangeDimensionalization(Vs,Ps,Es,NewParDef)
% NewVarDef is a matrix containing the relations of the new parameters to the old:
% All parameters are assumed to start consecutively after the first 2 functions

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

VsNew=zeros(size(Vs));
%Tnew = 0;

% Assume parameters start from the 4rd field in the Ps structure
defstart=4;
plen = size(NewParDef,2);

% Read parameters into an array
tempparms=[];
fn=fieldnames(Ps);
ind = defstart;
while length(tempparms)<plen
    tempparms = [tempparms; Ps.(fn{ind})(:)];
    ind=ind+1;
end;
pend=ind;

% Rescale the new vector using the old vector and transform matrix
newparms=ones(size(NewParDef,1),1);
for ind=1:length(tempparms);
    newparms=newparms.*(tempparms(ind).^NewParDef(:,ind));
end;


PsNew=Ps;
curloc=0;
% Put rescaled values into the correct parameters
for ind=defstart:pend-1
    %ind
    PsNew.(fn{ind})=newparms(curloc+(1:length(Ps.(fn{ind}))))';
    curloc=curloc+length(Ps.(fn{ind}));
end;

% Rescale variables
if(~isempty(Vs))
    if(length(newparms)<length(tempparms)+Ps.VarNum)
        warning('Could not rescale variables, Missing data in transform matrix.');
    else
        for ii=1:Ps.VarNum
            VsNew(:,ii,:)=Vs(:,ii,:)*newparms(length(tempparms)+ii);
        end;
    end;
end;

% Rescale space and time
if(length(newparms)>length(tempparms)+Ps.VarNum)
    PsNew.Lx = Ps.Lx * newparms(length(tempparms)+Ps.VarNum+1);
    PsNew.Ly = Ps.Ly * newparms(length(tempparms)+Ps.VarNum+1);
end;

scaling = newparms./[tempparms; ones(length(newparms)-length(tempparms),1)];

