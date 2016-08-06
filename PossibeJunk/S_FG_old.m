function Out=S_FG_old(Vs,Ps,Es,varargin)
% Spatial terms of Full Gilad model
% VsOut=S_FG(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the spatial terms of the model

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

if(Es.Jcob==0)
    if(Ps.Nx==1)
		len  = Ps.Ly;
        pnum = Ps.Ny;
	else	
		len  = Ps.Lx;
        pnum = Ps.Nx;
	end;

dx  = len/pnum;
xx  = ((1:pnum)*dx - len/2)';
top = (-xx.^2)/2;
lo2 = round(pnum/2);
%nrm = dx/(2*pi);   % for 2d
nrm = (dx/sqrt(2*pi))*(1+Ps.eta.*Vs(:,1)/5);    % for 1d

% gaus = exp(-xx.^2/2/(1+eta*b).^2)*w*dx;
%plot(exp(top/(1+7)^2))
%plot(1:length(xx),top)
OutB= zeros(pnum,1);
OutW= zeros(pnum,1);
%plot(circshift(exp(top./(1+Ps.eta*mean(Vs(:,1))).^2),lo2));
for ii=1:pnum
	% Calculate gaussian in 3 steps
	coef = Ps.eta*Vs(ii,1);
	gaus = circshift(exp(top./(1+coef).^2),ii+lo2);
	uptk = gaus.*Vs(:,2).*nrm*Vs(ii,1);
	%plot(uptk)
	%pause
	% update the biomass and water variables
	OutB(ii) = sum(uptk).*Ps.ni*(1-Vs(ii,1));   
	%size(OutB),size(OutW),size(uptk)
	OutW = OutW - uptk.*Ps.gamma;
end;
%plot([OutB-Ps.ni.*Vs(:,2).*Vs(:,1).*(1 + Ps.eta.*Vs(:,1)).^2.*(1 - Vs(:,1))])
%pause
%plot([OutW+Ps.gamma.*Vs(:,2).*(1 + Ps.eta.*Vs(:,1)).^2.*Vs(:,1)])
%pause
%subplot(1,2,1); plot(OutB);
%subplot(1,2,2); plot(OutW);
% Correct the "local parts" due to the form of L_SG file
% (not needed if using the "correct" kind)
%OutB=0; OutW=0;
OutB = OutB - Ps.ni.*Vs(:,2).*Vs(:,1).*(1 + Ps.eta.*Vs(:,1)).^2.*(1 - Vs(:,1));
OutW = OutW + Ps.gamma.*Vs(:,2).*(1 + Ps.eta.*Vs(:,1)).^2.*Vs(:,1);
%plot(OutW)
%pause
%plot(OutB)
%pause
if(length(Ps.SpaMat)<len)
    warning('No Spatial Matrix was defined');
    Ps.SpaMat = DervSM(2,Ps,Es);
end
% Update the diffusion parts
OutB = OutB + Ps.SpaMat*Vs(:,1).*Ps.Ds(1);
OutW = OutW + Ps.SpaMat*Vs(:,2).*Ps.Ds(2);
OutH = Ps.SpaMat*(Vs(:,3).^2).*Ps.Ds(3);
Out = [OutB OutW OutH];

else
warning('Jacobian not implemented');    
Out = 0;

end

