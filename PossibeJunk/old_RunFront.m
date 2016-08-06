function [points,VsOut]=RunFront(Vs1,Vs2,Ps,Es,IntFunc,TestFunc,varargin)
% Use an integrator and test function to run Front dynamics of two fileds (Vs1,Vs2).

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs1,Ps,Es,varargin{:});

% Check if output should be continously written out
Ps.Nx = Ps.Nx*2; Ps.Lx = Ps.Lx*2;

points=[];

VsFront = [Vs1;Vs2];
VsOut = RunToSS(VsFront,Ps,Es,IntFunc);
% Analysze the VsFront results
points = TestFunc(VsOut,Ps,Es);

end
