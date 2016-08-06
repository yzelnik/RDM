function Out=S_SG(Vs,Ps,Es)
% Spatial function of Simplfied Gilad model, or similar model with d2(H^2) term
% Out=S_SG(Vs,Ps,Es)
% Given the state variables (Vs) and parameters (Ps), calculate the spatial terms of the model
% The H variable is assumed to the be last one, unless specified with Ps.H2
% The array Ps.Ds is assumed to contain the diffusion rates of the different reactants 

if(~isfield(Ps,'NLD'))
    Ps.NLD=0;
end;
if(~isfield(Ps,'H2') || Ps.H2==0)
   Ps.H2=Ps.Vnum; % By defauly the H^2 field is assumed to be the last one
end;

if(isfield(Es,'SetupMode') && Es.SetupMode)
    % Pre caclculate spatial matrix, for future use
   Out = Ps;
   Out.Derv2Mat = DervSM(2,Ps,Es);
   if(Ps.NLD==0) % More direct way, just apply second derivative on H^2

   else          % Less direct, calculate SM including values of H
       disp('NLD !');
       
       %len=Ps.Nx*Ps.Ny;
       %Out.SpaMat = BuildDynamicSM(Out,Es);
       %dervh = Out.Derv2Mat*Vs(:,Ps.H2);
       %size(dervh)
       %Out.SpaMat((Ps.H2-1)*len+(1:len),(Ps.H2-1)*len+(1:len)) = dervh;
       Es.fmod=1;
       Es.SetupMode=0;
       Out.SpaMat=S_SG(Vs,Out,Es);
       Es.fmod=0;
       
   end;
else            % Normal run
   len=Ps.Nx*Ps.Ny; 
   
   if(~isfield(Ps,'Derv2Mat'))        % Caclculate spatial matrix if needed
        Ps.Derv2Mat = DervSM(2,Ps,Es);
   end;
   
   if(~isfield(Es,'fmod') || (Es.fmod==0))	% Model equations
        if(Ps.NLD==0)
            Out = [];
            for ind=1:Ps.Vnum
                if(ind==Ps.H2)
                    tmpout = Ps.Derv2Mat*(Vs(:,ind).^2).*Ps.Ds(ind);
                else
                    tmpout = Ps.Derv2Mat*Vs(:,ind).*Ps.Ds(ind);
                end;
            Out  = [Out tmpout];
            end;
        else % Use the jacobian as the SM
            Es.fmod=1;
            Out=S_SG(Vs,Ps,Es);
            Es.fmod=0;
        end;
   else			% Jacobian of equations
       
        Out = sparse(len*Ps.Vnum,len*Ps.Vnum);
        d2    = Ps.Derv2Mat;
        for ind=1:Ps.Vnum
            if(ind==Ps.H2)
                H     = Vs(:,ind);
                d1    = DervSM(1,Ps,Es);
        
                part1 = sparse(H*ones(1,len)).*d2;		% H0 * Derv2 (dH)
                part2 = sparse(diag(d2*H));			% dH * Derv2 (H0)
                part3 = 2*sparse((d1*H)*ones(1,len)).*d1;	% 2* Derv1(H0) * Derv1(dH)
                tmpout= 2* (part1 + part2 + part3);
            else
                tmpout= d2;
            end;
            Out(len*(ind-1)+(1:len),len*(ind-1)+(1:len)) = tmpout*Ps.Ds(ind);
        end;
        
   end;
end;            % End normal run

end

%%%%%%%%%%%%%%%%%% AUX FUNC %%%%%%%%%%%%%%%%%%

function sm = BuildDynamicSM(Ps,Es)
len=Ps.Nx*Ps.Ny;
for ii=1:Ps.Vnum	
    if(ii==Ps.H2)
        % Put in derivative sub-matrix in a block-diagonal fashion
        sm((ii-1)*len+(1:len),(ii-1)*len+(1:len)) = Ps.Derv2Mat*Ps.Ds(ii);
    end;
end;

end