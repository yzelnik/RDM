function MatOut=S_LD(~,Ps,Es)
% Spatial Matrix for a system with linear derivatives 
% MatOut=S_LD(Vs,Ps,Es)
% Ps.Nx (Ps.Ny),Ps.Lx (Ps.Ly) give the number of sites and system size for creating the matrix
% The array Ps.Ds is assumed to contain the coefficients for the different derivatives
% With the all the first derivatives first, then all second derivatives, and so on

if(isfield(Es,'SetupMode') && Es.SetupMode)
    % Pre caclculate spatial matrix, for future use
   MatOut = Ps;
   MatOut.SpaMat = CalcSM(Ps,Es);

else		% Normal run
   if(~isfield(Ps,'SpaMat'))    % Caclculate spatial matrix if needed
        Ps.SpaMat = CalcSM(Ps,Es);
   end;
   
   MatOut = Ps.SpaMat;
end

end

%------------ AUX funcs ----------------------------------------------------

function SM=CalcSM(Ps,Es)
   len=Ps.Nx*Ps.Ny;
   
   if(size(Ps.Ds,2)==1) Ps.Ds=Ps.Ds'; end; % Make sure Ds is in a row shape
   % find the highest derivative and organize the Ds vector of derivatives coefficients
   maxderv=ceil(size(Ps.Ds,2)/Ps.VarNum);
   if(size(Ps.Ds,2)<Ps.VarNum*maxderv) 
       Ps.Ds(1,Ps.VarNum*maxderv)=0;
   end;
%   Ds=reshape([Ps.Ds(:); zeros(maxderv*Ps.VarNum-length(Ps.Ds),1)],Ps.VarNum,maxderv);
   
   SM = sparse(len*Ps.VarNum,len*Ps.VarNum);	% Build the initial block-diagonal matrix
   for ii=1:maxderv
	if(sum(abs(Ps.Ds(:,ii))))			% If there are non zero coefficients for this derivative order
		base = DervSM(ii,Ps,Es);	% Get derivative spatial matrix
		for jj=1:Ps.VarNum
			loc = (jj-1)*len+(1:len); 			% Find location for this block-diagonal part
			SM(loc,loc) = SM(loc,loc) + base.*Ps.Ds(:,jj+(ii-1)*Ps.VarNum);	% Update this part
		end;
	end;
   end;
		
end
