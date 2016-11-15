function st=getode(Vs,Ps,Es,varargin)
% Get an ODE (non-spatial) solution, as defined by Ps.LocFunc
% st=getode(Vs,Ps,Es)

if(~mod(nargin,2)) error('No default extra-input exists for getode.'); end;

% Update online
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});
% Make sure Ps parameters are properly setup
[Vs,Ps,Es]=FillMissingPs(Vs,Ps,Es);

if(size(Vs,2)<Ps.VarNum)	% Replicate values, for the lazy amongst us
	tmp = repmat(Vs,1,Ps.VarNum);
	Vs = tmp(:,1:Ps.VarNum);
end;
    
% Find the steady-state
st=FindODESS(Vs,Ps,Es);

end