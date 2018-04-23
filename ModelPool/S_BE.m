function MatOut=S_BE(Vs,Ps,Es)
% Spatial function for Burger's equation
% MatOut=S_BE(Vs,Ps,Es)
% u_t = -d1 u*u_x + d2*u_xx

if(isfield(Es,'SetupMode') && Es.SetupMode)
   % Pre calculate spatial data, for future use
   MatOut = Ps;
   MatOut.SpaData = CalcSpaData(Ps,Es);
else		% Normal run
   
	len=Ps.Nx*Ps.Ny;
	if(~isfield(Ps,'SpaData'))    % calculate spatial data if needed
        Ps.SpaData = CalcSpaData(Ps,Es);
	end;
   
	if(~isfield(Es,'JacMode') || (Es.JacMode==0))	% Model equation
        MatOut = zeros(size(Vs,1),size(Vs,2)); 
        for ii=1:Ps.VarNum				
            MatOut(:,ii) = Ps.Ds(Ps.VarNum+ii)*Ps.SpaData{2}*Vs(:,ii) - Vs(:,ii).*(Ps.Ds(ii)*Ps.SpaData{1}*Vs(:,ii));
        end;
	else        % Jacobian of equation
        MatOut = zeros(size(Vs,1),size(Vs,2)); 
        for ii=1:Ps.VarNum				
            part1 = sparse(Vs(:,ii)*ones(1,len)).*Ps.SpaData{1};        % U0 * Derv1 (dU)
            part2 = sparse(diag(Ps.SpaData{1}*Vs(:,ii)));               % dU * Derv1 (U0)
            MatOut(len*(ii-1)+(1:len),len*(ii-1)+(1:len))= Ps.Ds(Ps.VarNum+ii)*Ps.SpaData{2}-Ps.Ds(ii)*(part1 + part2); 
        end;
   end;
end  % End normal run

end

%------------ AUX funcs ----------------------------------------------------

function SM=CalcSpaData(Ps,Es)	
   SM{1} = DervSM(1,Ps,Es); % Calc 1st derivative spatial matrix
   SM{2} = DervSM(2,Ps,Es); % Calc 2nd derivative spatial matrix
end

