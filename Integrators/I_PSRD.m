function VsOut=I_PSRD(Vs,Ps,Es,varargin)
% Integrator - Pseudo-Spectral method for Reaction-Diffusion
% Based on code by Aly-Khan Kassam (Solving reaction-diffusion equations 10 times faster)
% VsOut=I_PSRD(Vs,Ps,Es)
% The non-local part is assumed to be diffusion in periodic boundary conditions

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(isfield(Es,'SetupMode') && Es.SetupMode)
    % Setup variables for integration in the future    
    Ps=GetSpaData(Vs,Ps,Es);    % Run subfunction (see at the end)
    VsOut = Ps; % This is a "misuse" of the name, but we just want to return the "new" Ps struct
    
else        % Normal run
    if(~isfield(Ps,'SpaData') || isempty(Ps.SpaData) || ~isfield(Ps.SpaData,'h') || ~(Ps.SpaData.h==Es.Tstep) )  
        Ps=GetSpaData(Vs,Ps,Es);    % Run subfunction if needed (running this way is best avoided)
    end;
    
    % Initialization
    Vnum = Ps.Vnum;
    tmax = Es.Tdest;
    nmax = round(tmax/Ps.SpaData.h);
    tlen = Ps.Nx*Ps.Ny;
    
    for ii=1:Vnum
        Us{ii}=fftn(reshape(Vs(:,ii),Ps.Nx,Ps.Ny)); 
    end;
    %==================== TIME STEPPING LOOP =======================
    % The basic time-step loop
    for n = 1:nmax
        t=n*Ps.SpaData.h;
        
        
        for ii=1:Vnum   % Bring data back to real space and reshape
            Ts(:,ii)=reshape(ifftn(Us{ii}),tlen,1); 
        end;
        % Run local part
        TTs=Ps.LocFunc(Ts,Ps,Es);
        % Run first spatial part
        for ii=1:Vnum  
            Nvs{ii}=fftn(reshape(TTs(:,ii),Ps.Nx,Ps.Ny));		%Nonlinear evaluation. g(u,*)
            As{ii}=Ps.SpaData.EE1s{ii}.*Us{ii} + Ps.SpaData.Qs{ii}.*Nvs{ii};          %Coefficient 'a' in ETDRK formula 	
            Ts(:,ii)=reshape(ifftn(As{ii}),tlen,1);             % Reshape before running local part
        end;
        % Run local part
        TTs=Ps.LocFunc(Ts,Ps,Es);
        % Run second spatial part
        for ii=1:Vnum
            Nas{ii}=fftn(reshape(TTs(:,ii),Ps.Nx,Ps.Ny));		%Nonlinear evaluation. g(a,*)
            Bs{ii}=Ps.SpaData.EE1s{ii}.*Us{ii} + Ps.SpaData.Qs{ii}.*Nvs{ii};          %Coefficient 'b' in ETDRK formula
            Ts(:,ii)=reshape(ifftn(Bs{ii}),tlen,1);             % Reshape before running local part
        end;
        % Run local part
        TTs=Ps.LocFunc(Ts,Ps,Es);
        % Run third spatial part
        for ii=1:Vnum
            Nbs{ii}=fftn(reshape(TTs(:,ii),Ps.Nx,Ps.Ny));			%Nonlinear evaluation. g(b,*)
            Cs{ii}=Ps.SpaData.EE1s{ii}.*As{ii} + Ps.SpaData.Qs{ii}.*(2*Nbs{ii}-Nvs{ii});	%Coefficient 'c' in ETDRK formula
            Ts(:,ii)=reshape(ifftn(Cs{ii}),tlen,1);                 % Reshape before running local part
        end;
        % Run local part
        TTs=Ps.LocFunc(Ts,Ps,Es);
        % Run final spatial part
        for ii=1:Vnum
            Ncs{ii}=fftn(reshape(TTs(:,ii),Ps.Nx,Ps.Ny));			%Nonlinear evaluation. g(c,*)
            Us{ii}=Ps.SpaData.E1s{ii}.*Us{ii} + Nvs{ii}.*Ps.SpaData.f1s{ii} + (Nas{ii}+Nbs{ii}).*Ps.SpaData.f2s{ii} + Ncs{ii}.*Ps.SpaData.f3s{ii};	%update
            Us{ii}(Ps.SpaData.Dealias_ind) = 0;					% High frequency removal --- de-aliasing
        end;
    end

    % Bring back to real space
    for ii=1:Vnum
        VsOut(:,ii)=reshape(real(ifftn(Us{ii})),tlen,1); 
    end;

    if(~isfield(Es,'SkipWarning') || Es.SkipWarning==0)    % Check for problematic values
        if any(isnan(reshape(VsOut,1,numel(VsOut))))
            warning('Matrix contains NaN')
        end;
    end;
end;

end

%%%%%%%%%%%%%%%%%  AUX function to prep things before integration %%%%%%%%%%%%%%%%%  
function Ps=GetSpaData(Vs,Ps,Es)
    %disp('running AUX GetSpaData function');
    Ps.SpaData.h    = Es.Tstep;
    % Dimensions
    resx = Ps.Lx/Ps.Nx;
    resy = Ps.Ly/Ps.Ny;
    if(Ps.Ly==0)  
        resy  = rex; 
        Ps.Ly = resy; 
    end;

    %========================= CUSTOM SET UP =========================
    x=resx*(1:Ps.Nx)'; y=resy*(1:Ps.Ny)'; [X,Y]=ndgrid(x,y);	%{[X,Y,Z]=ndgrid(x,x,x)}
    kx=[0:Ps.Nx/2-1 0 -Ps.Nx/2+1:-1]'/(Ps.Lx/(2*pi)); 		%wave numbers
    ky=[0:Ps.Ny/2-1 0 -Ps.Ny/2+1:-1]'/(Ps.Ly/(2*pi));
    [xi,eta]=ndgrid(kx,ky); 			%2D wave numbers. {[xi,eta,zeta]=ndgrid(k,k,k)}

    if(length(Ps.Ds)==Ps.Vnum)          % Assume Ds contains only 2nd derivatives if the size of Ds is appropiate
        justdiffusion=1;
    else
        justdiffusion=0;
    end;
    Ds   = [Ps.Ds(:) ;zeros(16,1)];
    for ii=1:Ps.Vnum
        Us{ii}=fftn(reshape(Vs(:,ii),Ps.Nx,Ps.Ny));
        if(justdiffusion==1)
            Ls{ii}=-Ds(ii)*(eta.^2+xi.^2);
        else
            Ls{ii}=Ds(ii)*i*xi;                                 %2D Grad(x)
            Ls{ii}=Ls{ii} - Ds(Ps.Vnum+ii)*(eta.^2+xi.^2);         %2D Laplacian. {-D*(eta.^2+xi.^2+zeta.^2)}
            Ls{ii}=Ls{ii} - Ds(Ps.Vnum*2+ii)*i*xi.*(eta.^2+xi.^2);	%2D Laplacian*Grad(x)
            Ls{ii}=Ls{ii} + Ds(Ps.Vnum*3+ii)*(eta.^2+xi.^2).^2;	%2D Laplacian^2
        end;
    end;

    Frx=logical(zeros(Ps.Nx,1));				%High frequencies for de-aliasing
    Frx([Ps.Nx/2+1-round(Ps.Nx/6) : Ps.Nx/2+round(Ps.Nx/6)])=1;
    Fry=logical(zeros(Ps.Ny,1));				%High frequencies for de-aliasing
    Fry([Ps.Ny/2+1-round(Ps.Ny/6) : Ps.Ny/2+round(Ps.Ny/6)])=1;
    [alxi,aleta]=ndgrid(Frx,Fry);				%{[alxi,aleta,alzeta]=ndgrid(Fr,Fr,Fr)}
    Ps.SpaData.Dealias_ind=alxi | aleta; 					%de-aliasing index. {alxi | aleta | alzeta}

    %=============== PRECOMPUTING ETDRK4 COEFFS =====================
    M=16; 							% no. of points for complex mean
    r=exp(i*pi*((1:M)-0.5)/M); 				% roots of unity

    % For each variable (with seperate diffusion rate) pre-compute coefficients:
    for ii=1:Ps.Vnum
        Ps.SpaData.E1s{ii}=exp(Ps.SpaData.h*Ls{ii}); 
        Ps.SpaData.EE1s{ii}=exp(Ps.SpaData.h*Ls{ii}/2);
        Ls{ii} =Ls{ii}(:); 
        LR=Ps.SpaData.h*Ls{ii}(:,ones(M,1))+r(ones(Ps.Nx*Ps.Ny,1),:);			%{r(ones(N^3,1),:)}
        Ps.SpaData.Qs{ii} = Ps.SpaData.h*real(mean( (exp(LR/2)-1)./LR					,2));
        Ps.SpaData.f1s{ii}=Ps.SpaData.h*real(mean( (-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3 	,2));
        Ps.SpaData.f2s{ii}=Ps.SpaData.h*real(mean( (4+2*LR+exp(LR).*(-4+2*LR))./LR.^3 		,2));
        Ps.SpaData.f3s{ii}=Ps.SpaData.h*real(mean( (-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3 	,2));
        Ps.SpaData.f1s{ii}=reshape(Ps.SpaData.f1s{ii},Ps.Nx,Ps.Ny); 
        Ps.SpaData.f2s{ii}=reshape(Ps.SpaData.f2s{ii},Ps.Nx,Ps.Ny); 
        Ps.SpaData.f3s{ii}=reshape(Ps.SpaData.f3s{ii},Ps.Nx,Ps.Ny);
        Ps.SpaData.Qs{ii}=reshape(Ps.SpaData.Qs{ii},Ps.Nx,Ps.Ny); 		%{reshape(*,N,N,N)}
    end;
end
