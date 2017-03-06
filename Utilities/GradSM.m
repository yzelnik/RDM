function gradient=GradSM(order,Ps,Es,varargin)
% Build a gradient operator using first derivative spatial matrices
% Matt=GradSM(order,Ps,Es)

% Update online if necessary
if(nargin>3) [~,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:}); end;

% Is this a 1D system?
if((Ps.Nx==1) || (Ps.Ny==1))
    gradient={DervSM(1,Ps,Es)};
  else  % 2D system
    gradient={DervSM([1 1],Ps,Es),DervSM([1 2],Ps,Es)};  
end;

end