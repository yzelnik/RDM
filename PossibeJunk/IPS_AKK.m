function [VsOut,uu]=IPS_AKK(Vs,Ps,Es,varargin)
% Integrator - Pseudo-Spectral, based on code by Aly-Khan Kassam (Solving reaction-diffusion equations 10 times faster)
% [VsOut,uu]=IPS_AKK(Vs,Ps,Es)
% The non-local part is assumed to be a diffusion in periodic boundary conditions

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

%========================= CUSTOM SET UP =========================
% Initialization
Vnum = Ps.Vnum;
Ds   = Ps.Ds;
h    = Es.Tstep;
tmax = Es.Tdest;

% Dimensions
Lx   = Ps.Lx;
Ly   = Ps.Ly;
Nx   = Ps.Nx;
Ny   = Ps.Ny;
resx = Lx/Nx;
resy = Ly/Ny;
if(Ly==0)  resy=rex; Ly=resy; end;


x=resx*(1:Nx)'; y=resy*(1:Ny)'; [X,Y]=ndgrid(x,y);	%{[X,Y,Z]=ndgrid(x,x,x)}
kx=[0:Nx/2-1 0 -Nx/2+1:-1]'/(Lx/(2*pi)); 		%wave numbers
ky=[0:Ny/2-1 0 -Ny/2+1:-1]'/(Ly/(2*pi));
[xi,eta]=ndgrid(kx,ky); 				%2D wave numbers. {[xi,eta,zeta]=ndgrid(k,k,k)}


for ii=1:Vnum
	Us{ii}=fftn(reshape(Vs(:,ii),Nx,Ny)); 
	Ls{ii}=-Ds(ii)*(eta.^2+xi.^2);			%2D Laplacian. {-D*(eta.^2+xi.^2+zeta.^2)}
end;
			
Frx=logical(zeros(Nx,1));				%High frequencies for de-aliasing
Frx([Nx/2+1-round(Nx/6) : Nx/2+round(Nx/6)])=1;
Fry=logical(zeros(Ny,1));				%High frequencies for de-aliasing
Fry([Ny/2+1-round(Ny/6) : Ny/2+round(Ny/6)])=1;
[alxi,aleta]=ndgrid(Frx,Fry);				%{[alxi,aleta,alzeta]=ndgrid(Fr,Fr,Fr)}
ind=alxi | aleta; 					%de-aliasing index. {alxi | aleta | alzeta}

%=============== PRECOMPUTING ETDRK4 COEFFS =====================
M=16; 							% no. of points for complex mean
r=exp(i*pi*((1:M)-0.5)/M); 				% roots of unity

% For each variable (with seperate diffusion rate) pre-compute coefficients:
for ii=1:Vnum
	E1s{ii}=exp(h*Ls{ii}); 
	EE1s{ii}=exp(h*Ls{ii}/2);
	Ls{ii}=Ls{ii}(:); 
	LRs{ii}=h*Ls{ii}(:,ones(M,1))+r(ones(Nx*Ny,1),:);			%{r(ones(N^3,1),:)}
	Qs{ii}= h*real(mean( (exp(LRs{ii}/2)-1)./LRs{ii}					,2));
	f1s{ii}=h*real(mean( (-4-LRs{ii}+exp(LRs{ii}).*(4-3*LRs{ii}+LRs{ii}.^2))./LRs{ii}.^3 	,2));
	f2s{ii}=h*real(mean( (4+2*LRs{ii}+exp(LRs{ii}).*(-4+2*LRs{ii}))./LRs{ii}.^3 		,2));
	f3s{ii}=h*real(mean( (-4-3*LRs{ii}-LRs{ii}.^2+exp(LRs{ii}).*(4-LRs{ii}))./LRs{ii}.^3 	,2));
	f1s{ii}=reshape(f1s{ii},Nx,Ny); 
	f2s{ii}=reshape(f2s{ii},Nx,Ny); 
	f3s{ii}=reshape(f3s{ii},Nx,Ny);
	Ls{ii}=reshape(Ls{ii},Nx,Ny); 
	Qs{ii}=reshape(Qs{ii},Nx,Ny); 		%{reshape(*,N,N,N)}
end;
clear LRs

%==================== TIME STEPPING LOOP =======================
nmax=round(tmax/h);
% The basic time-step loop
for n = 1:nmax
	t=n*h;
	% Bring it back to real space
	for ii=1:Vnum
		Ts(:,ii)=reshape(ifftn(Us{ii}),Nx*Ny,1); 
	end;
	% Run local part
	TTs=Ps.LocFunc(Ts,Ps,Es);
	% Run first spatial part
	for ii=1:Vnum  
		Nvs{ii}=fftn(reshape(TTs(:,ii),Nx,Ny));		%Nonlinear evaluation. g(u,*)
		As{ii}=EE1s{ii}.*Us{ii} + Qs{ii}.*Nvs{ii};	%Coefficient ’a’ in ETDRK formula 	
		Ts(:,ii)=reshape(ifftn(As{ii}),Nx*Ny,1); 
	end;
	% Run local part
	TTs=Ps.LocFunc(Ts,Ps,Es);
	% Run second spatial part
	for ii=1:Vnum
		Nas{ii}=fftn(reshape(TTs(:,ii),Nx,Ny));		%Nonlinear evaluation. g(a,*)
		Bs{ii}=EE1s{ii}.*Us{ii} + Qs{ii}.*Nas{ii};	%Coefficient ’b’ in ETDRK formula
		Ts(:,ii)=reshape(ifftn(Bs{ii}),Nx*Ny,1); 
	end;
	% Run local part
	TTs=Ps.LocFunc(Ts,Ps,Es);
	% Run third spatial part
	for ii=1:Vnum
		Nbs{ii}=fftn(reshape(TTs(:,ii),Nx,Ny));			%Nonlinear evaluation. g(b,*)
		Cs{ii}=EE1s{ii}.*As{ii} + Qs{ii}.*(2*Nbs{ii}-Nvs{ii});	%Coefficient ’c’ in ETDRK formula
		Ts(:,ii)=reshape(ifftn(Cs{ii}),Nx*Ny,1); 
	end;
	% Run local part
	TTs=Ps.LocFunc(Ts,Ps,Es);
	% Run final spatial part
	for ii=1:Vnum
		Ncs{ii}=fftn(reshape(TTs(:,ii),Nx,Ny));			%Nonlinear evaluation. g(c,*)
		Us{ii}=E1s{ii}.*Us{ii} + Nvs{ii}.*f1s{ii} + (Nas{ii}+Nbs{ii}).*f2s{ii} + Ncs{ii}.*f3s{ii};	%update
		Us{ii}(ind) = 0;					% High frequency removal --- de-aliasing
	end;

end
% Bring back to real space
for ii=1:Vnum
	VsOut(:,ii)=reshape(real(ifftn(Us{ii})),Nx*Ny,1); 
end;

if any(isnan(reshape(VsOut,1,numel(VsOut))))
	if(~isfield(Es,'SkipWarning') || Es.SkipWarning)
		warning('Matrix contains NaN')
	end
end

end
