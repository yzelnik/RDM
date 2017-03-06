function VsOut=I_FDCN(Vs,Ps,Es,varargin)
% Integrator with Finite-Difference Crank-Nicolson Method (semi-implicit)
% VsOut=I_FDCN(Vs,Ps,Es)
% Ps.LocFunc and Ps.SpaFunc are the functions that detail the local and non-local parts of the system

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(isfield(Es,'SetupMode') && Es.SetupMode)
    % Setup variables for integration in the future    
    Ps=GetSpaData(Vs,Ps,Es);    % Run subfunction (see at the end)
    VsOut = Ps; % This is a "misuse" of the name, but we just want to return the "new" Ps struct
    
else            % Normal run
    % Setup the spatial matrix and auxiliary flags if not already done
    if(~isfield(Ps,'SpaMat') || isempty(Ps.SpaMat))
        [Vs,Ps,Es]=SetupSpatialData(Vs,Ps,Es);
    end;
    
    if(~isfield(Ps,'SpaData') || isempty(Ps.SpaData) || ~isfield(Ps.SpaData,'TsSize') || ~(Ps.SpaData.TsSize==Es.TsSize) )  
        Ps=GetSpaData(Vs,Ps,Es);    % Run subfunction if needed (running this way is best avoided)
    end;
    
    if(~Es.SmUse)
        error('Implicit integration is not possible without a matrix strucutre for the spatial part of the PDE.');
    end;
    
    NonNeg = 0;		% If we know variables are positive, make sure they remain so	
    if((isfield(Es,'NonNeg')) & (Es.NonNeg))
        NonNeg = 1;
    end;
    
    %tempmat = speye(size(Ps.SpaMat)) - Es.TsSize*Ps.SpaMat;
    syslen = Ps.Nx * Ps.Ny;
    totsteps = ceil(Es.TimeDst/Es.TsSize);
    % Go through each time step
    for ii=1:totsteps 
        % Explicit step for local (non-linear) part
        VsLocChange = reshape( Es.TsSize*Ps.LocFunc(Vs,Ps,Es) ,syslen*Ps.VarNum,1);	
        % Implicit step for spatial (potentially linear) part
        VsTemp = Ps.SpaData.InvIminusSM * ( Ps.SpaData.IplusSM*reshape(Vs,syslen*Ps.VarNum,1)  + VsLocChange );    
        Vs     = reshape(VsTemp,syslen,Ps.VarNum);
        
        if Es.SmUpdate     % Use this if the spatial matrix needs to be updated online
            Ps=GetSpaData(Vs,Ps,Es);    % update the necessary spatial information for this integration
        end;
        
        if(NonNeg)  % consider deleting
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