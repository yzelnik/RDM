function [points,VsFrontTotal,VsFrontOutTotal]=ContParFront(Vs,VsHom,Ps,Es,IntFunc,TestFunc,varargin)
% Use an integrator and test function for continuation of a paramterer
% For each value of continuation parameter, uniform and nonuniform solutions
% are obtained, and the dynamics of a front between the two is simulated.
% Vs = nonuniform solution; VsHom = uniform (homogeneous) solution;

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

% Check if output should be continously written out
WriteFlag = 0;
if(isfield(Es,'BFout') & Es.BFout)
	WriteFlag = 1;
end;
parname = Es.BFpar;
parrange = Es.BFrange;

PsHom = Ps; PsFront = Ps;
PsHom.Ds = ones(size(PsHom.Ds)).*Es.STsmall;
PsFront.Nx = Ps.Nx*2; PsFront.Lx = Ps.Lx*2;

EsHom = Es; EsFront = Es;
EsHom.STsmall = 1e-4;
EsFront.STsmall = 1e-4;

points=[];
VsFrontTotal = zeros(PsFront.Nx,Ps.Vnum,length(parrange));
VsFrontOutTotal = zeros(PsFront.Nx,Ps.Vnum,length(parrange));

for ii=1:length(parrange)
	% Update paramater
	Ps.(parname) = parrange(ii);
	PsHom.(parname) = parrange(ii);
	PsFront.(parname) = parrange(ii);
    disp([Es.BFpar '=' num2str(parrange(ii))]);

	% Run both nonuniform and uniform systems to SS
	VsHomOut = RunToSS(VsHom,PsHom,EsHom,IntFunc);
	Vs = Vs + rand(size(Vs))*Es.STsmall;
	VsOut = RunToSS(Vs,Ps,Es,IntFunc);
	% Analysze the nonuniform results
	ResVs = TestFunc(VsOut,Ps,Es);

	% prepare and run front between uniform and nonuniform
	VsFront = [VsOut;VsHomOut];
	VsFrontOut = RunToSS(VsFront,PsFront,EsFront,IntFunc);
	% Analysze the VsFront results
	ResVsFront = TestFunc(VsFrontOut,PsFront,EsFront);

	% Add results to the overall results
	points = [points; Ps.(parname) ResVs(:)' ResVsFront(:)'];
	if(WriteFlag)
		dlmwrite(Es.BFout,points);
	end;
	VsFrontTotal(:,:,ii) = VsFront;
	VsFrontOutTotal(:,:,ii) = VsFrontOut;

	% Use the last result to get faster convergence
	Vs = VsOut;
	VsHom = VsHomOut;

end;
