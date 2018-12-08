function Out=S_FG(Vs,Ps,Es)
% Spatial Matrix of Full Gilad model
% VsOut=S_FG(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the spatial terms of the model

if(Ps.Nx==1)	% This handles only for 1D, should be changed...
   len  = Ps.Ly;
   pnum = Ps.Ny;
elseif(Ps.Ny==1)
   len  = Ps.Lx;
   pnum = Ps.Nx;
else
   error('Only 1D supported');
end;

if(isfield(Es,'SetupMode') && Es.SetupMode)
    % Pre caclculate spatial matrix and gaussian, for future use
   Out = Ps;
   Out.Derv2Mat = DervSM(2,Ps,Es);
   Out.GausNM = CalcGausNM(Ps,Es);

else		% Normal run
    % Check for missing fields
   if((~isfield(Ps,'GausNM')) || (~isfield(Ps,'Derv2Mat'))  || isempty(Ps.Derv2Mat) || isempty(Ps.GausNM))              
        Ps.Derv2Mat = DervSM(2,Ps,Es);
        Ps.GausNM = CalcGausNM(Ps,Es);
   end;
   B = Vs(:,1);
   W = Vs(:,2);
   
   dx  = len/pnum; 
   nrm = dx/(2*pi*Ps.S0^2);   % for 2d
   %nrm = (dx/sqrt(2*pi*Ps.S0^2))*repmat(1+Ps.eta.*Vs(:,1),1,pnum);    % for 1d

   if(~isfield(Es,'JacMode') || (Es.JacMode==0))  % Model equations
	
        % Calculate a matrix of "connections" between each B and W in space
        coef = repmat((Ps.eta.*B)',pnum,1); 	% calculate the coefficient (eta*B) in a matrix form
        gaus = exp(Ps.GausNM./(1+coef).^2);         % insert coef & pre-calcualted Gaussian Nominator Matrix into gaussian
        uptk = gaus.*repmat(W,1,pnum).*repmat(B',pnum,1).*nrm;	% Calculate the uptake matrix using the gaussian

        OutB = sum(uptk,1)'.*Ps.lambda.*(1-B);	% Sum columns, and mix in the carrying capacity, to get effect on B
        OutW = -sum(uptk,2).*Ps.gamma;              % Sum rows to get effect on W

        % Correct the "local parts" due to the form of L_SG file (not needed if using the "correct" kind)
        SimpB= Ps.lambda.*W.*B.*(1 + Ps.eta.*B).^2.*(1 - B);
        SimpW= Ps.gamma.*W.*(1 + Ps.eta.*B).^2.*B;
        OutB = OutB - SimpB;
        OutW = OutW + SimpW;

        % Update the diffusion parts
        OutB = OutB + Ps.Derv2Mat*B.*Ps.Ds(1);
        OutW = OutW + Ps.Derv2Mat*W.*Ps.Ds(2);
        OutH = Ps.Derv2Mat*(Vs(:,3).^2).*Ps.Ds(3);
	
        Out = [OutB OutW OutH];		% Add them all up
	
   else			% Jacobian of equations
        warning('Jacobian not implemented');    
        Out = 0;
   end;   

end		% End normal run

end

%------------ AUX funcs ----------------------------------------------------

function SM=CalcGausNM(Ps,Es)	% Gaussian Nominator Matrix

if(Ps.Nx==1)	% This handles only for 1D, should be changed...
   len  = Ps.Ly;
   pnum = Ps.Ny;
elseif(Ps.Ny==1)
   len  = Ps.Lx;
   pnum = Ps.Nx;
end;

dx  = len/pnum; 
xx  = ((1:pnum)*dx - len/2)'; 
top = (-xx.^2)/(2*Ps.S0^2); 
lo2 = round(pnum/2);
for ii=1:pnum 
	SM(:,ii)=circshift(top,ii+lo2); 
end; 
		
end

