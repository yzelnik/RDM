function VsOut=I_FDCN(Vs,Ps,Es,varargin)
% Integrator with Finite-Difference Crank-Nicolson Method (semi-implicit)
% VsOut=I_FDCN(Vs,Ps,Es)
% Ps.LocFunc and Ps.SpaFunc are the functions that detail the local and non-local parts of the system

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Es=InsertDefaultValues(Es,'NonNeg',0,'SetupMode',0);

if(Es.SetupMode)
    % Setup variables for integration in the future    
    Ps=GetSpaData(Vs,Ps,Es);    % Run subfunction (see at the end)
    VsOut = Ps; % This is a "misuse" of the name, but we just want to return the "new" Ps struct
    
else            % Normal run
    % Setup the spatial matrix and auxiliary flags if not already done
    if(~isfield(Es,'SpaMatUse') || Es.SpaMatUse<0)
        [Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);
    end;
    
    if(~isfield(Ps,'SpaData') || isempty(Ps.SpaData) || ~isfield(Ps.SpaData,'TsSize') || ~(Ps.SpaData.TsSize==Es.TsSize) )  
        Ps=GetSpaData(Vs,Ps,Es);    % Run subfunction if needed (running this way is best avoided)
    end;
    
    if(~Es.SpaMatUse)
        error('Implicit integration is not possible without a matrix strucutre for the spatial part of the PDE.');
    end;
    
    %tempmat = speye(size(Ps.SpaMat)) - Es.TsSize*Ps.SpaMat;
    syslen = Ps.Nx * Ps.Ny;
    totsteps = ceil(Es.TimeDst/Es.TsSize);
    % Go through each time step
    for ii=1:totsteps 
        % Explicit step for local (probably non-linear) part
        VsLocChange = reshape( Es.TsSize*Ps.LocFunc(Vs,Ps,Es) ,syslen*Ps.VarNum,1);	
        % Implicit step for spatial (and linear) part
        VsTemp = Ps.SpaData.InvIminusSM * ( Ps.SpaData.IplusSM*reshape(Vs,syslen*Ps.VarNum,1)  + VsLocChange );    
        Vs     = reshape(VsTemp,syslen,Ps.VarNum);
        
        if(Es.NonNeg)  % make sure values are not negative, if relevant
            Vs = max(0,Vs);
        end;
    end; 
    VsOut = Vs;
end;        % End of normal run

end


%%%%%%%%%%%%%%%%%  AUX function to prep things before integration %%%%%%%%%%%%%%%%%  
function Ps=GetSpaData(~,Ps,Es)
    %disp('running AUX GetSpaData function');
    Ps.SpaData.IplusSM = speye(size(Ps.SpaMat)) + Es.TsSize/2*Ps.SpaMat;
    Ps.SpaData.InvIminusSM = inv(speye(size(Ps.SpaMat)) - Es.TsSize/2*Ps.SpaMat);
    Ps.SpaData.TsSize = Es.TsSize;
end