function MatOut=S_LL(Vs,Ps,Es)
% Spatial Matrix of Lejeune-Lefever model
% MatOut=S_LL(Vs,Ps,Es)
% Ps.Nx (Ps.Ny),Ps.Lx (Ps.Ly) give the number of sites and system size for creating the matrix

if(isfield(Es,'SetupMode') && Es.SetupMode)
    % Pre caclculate spatial matrix, for future use
   MatOut = Ps;
   MatOut.SpaData = CalcSpaData(Ps,Es);

else		% Normal run
   len=Ps.Nx*Ps.Ny;
   B = Vs(:,1,1);
   if(~isfield(Ps,'SpaData'))    % calculate spatial matrices if needed
        Ps.SpaData = CalcSpaData(Ps,Es);
   end;
   
   dsm2 = Ps.SpaData{1};   
   dsm4 = Ps.SpaData{2};
   
   if(~isfield(Es,'JacMode') || (Es.JacMode==0))	% Model equation
        MatOut = 0.5*(Ps.L^2 - B).*(dsm2*B) - 0.125*B.*(dsm4*B) ;
   else             % Jacobian of equation
	
        MatOut = 0.5*dsm2*Ps.L^2;                   % Start with the normal diffusion part
	
        part1 = sparse(B*ones(1,len)).*dsm2;        % B0 * Derv2 (dB)
        part2 = sparse(diag(dsm2*B));               % dB * Derv2 (B0)
        MatOut = MatOut - 0.5*(part1 + part2);      % Add the non-linear diffusion part
	
        part1 = sparse(B*ones(1,len)).*dsm4;        % B0 * Derv4 (dB)
        part2 = sparse(diag(dsm4*B));               % dB * Derv4 (B0)
        MatOut = MatOut - 0.125*(part1 + part2);	% Add the non-linear 4th derivative part
   end
end		% End normal run

end

%------------ AUX funcs ----------------------------------------------------

function SM=CalcSpaData(Ps,Es)	
   SM{1} = DervSM(2,Ps,Es); % Calc 2nd derivative spatial matrix
   SM{2} = DervSM(4,Ps,Es); % Calc 4th derivative spatial matrix
end

