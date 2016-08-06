function PsOut=ConvPars_DC(Vs,Ps,Es,varargin)
% Convert parameters for DC model
% Use rescaling to convert dimensional to non-dimensional parameters
% and create a new Ps structure
% Ps_DCn=ConvPars_DC([],Ps_DC,[]);

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

PsOut.LocFunc 	= @L_DCn;
PsOut.SpaFunc   = Ps.SpaFunc;
PsOut.p			= (Ps.LambdaB.*Ps.P)./(Ps.N.*Ps.MB);
PsOut.nu		= Ps.N./Ps.MB;
PsOut.eta		= Ps.EE.*Ps.KB;
PsOut.phib		= Ps.PhiB./Ps.MB;
PsOut.k			= (Ps.KC.^Ps.nn)./(Ps.KB.^Ps.mm);
PsOut.bZero		= Ps.BZero./Ps.KB;
PsOut.nn		= Ps.nn;
PsOut.mm		= Ps.mm;
PsOut.lambdacw	= Ps.LambdaCW./Ps.LambdaB;
PsOut.lambdach 	= Ps.LambdaCH./Ps.LambdaB;
PsOut.mu 		= Ps.MC./Ps.MB;
PsOut.phic		= Ps.PhiC.*Ps.KB./Ps.MB;
PsOut.alpha		= Ps.A/Ps.MB;
PsOut.qc		= Ps.QC./Ps.KC;
PsOut.fc		= Ps.fC;
PsOut.qb		= Ps.QB./Ps.KB;
PsOut.fb		= Ps.fB;
PsOut.r			= Ps.R;
PsOut.gammab	= Ps.GammaB.*Ps.KB./Ps.MB;
PsOut.gammacw	= Ps.GammaCW.*Ps.KC./Ps.MB;
PsOut.gammach 	= Ps.GammaCH.*Ps.KC./Ps.MB;
PsOut.Vnum      = Ps.Vnum;
PsOut.Ds 		= Ps.Ds./Ps.Ds(1);
PsOut.Lx		= Ps.Lx*sqrt(Ps.MB/Ps.Ds(1));
PsOut.Ly		= Ps.Ly*sqrt(Ps.MB/Ps.Ds(1));
PsOut.Nx        = Ps.Nx;
PsOut.Ny        = Ps.Ny;

end
