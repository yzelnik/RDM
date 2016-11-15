function VsOut=L_CVS(Vs,Ps,Es,varargin)
% Crust-Vegetation-Space (CVS)
% VsOut=L_CVS(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the local terms of the model
% Variables are: V=vegetation,B=crust

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

% Initialization

V=Vs(:,1); 
B=Vs(:,2); 

if (Ps.P >= Ps.PvMin)
    Ps.AlphaV = Ps.AlphaVMax*(1-exp(-(Ps.P-Ps.PvMin)/Ps.Cv));
else
    Ps.AlphaV = 0;
end

if (Ps.P >= Ps.PbMin)
    Ps.AlphaB = Ps.AlphaBMax*(1-exp(-(Ps.P-Ps.PbMin)/Ps.Cb));
else
    Ps.AlphaB = 0;
end

if(Es.JacMode==0)
	% Model equations
	g = 0.5.*(tanh(Ps.d.*(Ps.Vc-V))+1);

	dV = Ps.AlphaV.*(V+Ps.EtaV).*(1-V-B) - Ps.EpsV.*Ps.Dp.*g.*(1-V-B).*V - Ps.Gamma.*Ps.Dp^(2/3).*V - Ps.MuV.*V - Ps.PhiV.*V.*B;

	dB = Ps.AlphaB.*(B+Ps.EtaB).*(1-V-B) - Ps.EpsB.*Ps.Dp.*g.*(1-V-B).*B - Ps.MuB.*B - Ps.PhiB.*V.*B;

	VsOut = [dV,dB];

else
% Jacobian of equations
VdV = 0.*V;
VdB = 0.*V;
BdV = 0.*V;
BdB = 0.*V;

disp('not implemented yet');

% written in a large sparse matrix format 
VsOut = sparse([diag(VdV) diag(VdB); diag(BdV) diag(BdB)]);
end;


