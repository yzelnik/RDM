function VsNew=TransformVar(Vs,Ps,Es,NewVarDef,varargin)
% Transform a state in a given variable definition to a new definition
% VsNew=TransformVar(Vs,Ps,Es,NewVarDef)
% NewVarDef is a vector containing the definition of the new variables:
% val = i*Ps.Vnum+j means an i'th order derivative of the j'th original variable

% Update online if necessary
[Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:});

VsNew = zeros(size(Vs,1),length(NewVarDef),size(Vs,3));

if(min(NewVarDef)>0)  % Use only derivative transformations. 
    maxorder = ceil(max(NewVarDef)/Ps.Vnum)-1;
    dermats{1} = eye(size(Vs,1));  % dummy matrix - 0'th order derivative
    for ii=1:maxorder % matrix for ii'th order derivatives
        %size(dermats)
        dermats{ii+1} = DervSM(ii,Ps,Es);
    end;
    
    for ii=1:size(Vs,3)
        for varind = 1:length(NewVarDef)
            orgvar = mod(NewVarDef(varind)-1,Ps.Vnum)+1;
            orderd  = (NewVarDef(varind)-orgvar)/Ps.Vnum;
            VsNew(:,varind,ii) = dermats{orderd+1}*Vs(:,orgvar,ii);
            %disp([NewVarDef(varind) orgvar order]);
        end;
    end;
else
    warning('Non-derivative transformations Not Supported');
end;


end

