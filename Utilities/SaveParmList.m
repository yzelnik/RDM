function [Vs,Ps,Es]=SaveParmList(Vs,Ps,Es,vallist,namelist)
% Save a list of parameters into Ps/Es
% [Vs,Ps,Es]=SaveParmList(Vs,Ps,Es,vallist,namelist)
% vallist is a vector of parameter values (otherwise Es.ParmList is used)
% namelist is a cell array of parameter names (otherwise Es.BFpar is used)
% If Es.ParInCell exists and is non-zero, then it points to values as Es.CellData
% A proper text name is assumed to be a Ps variable (unless it is "Vs")
% A numerical value acts as a pointer in the Ps structure (offset by +3)
% Anything else is directly evaluated using eval (not as safe) such as: Ds(2)
if(nargin<4) vallist = Es.ParmList; end;
if(nargin<5) namelist = Es.BFpar; end;

for jj=1:length(namelist)
	if(isfield(Es,'ParInCell') && Es.ParInCell(jj))
        tmpval=Es.CellData{Es.ParInCell(jj)}{vallist(jj)};
    else
        tmpval=vallist(jj);
	end;
	if(isnumeric(namelist{jj})) % Allow access to model-parameters by index
        tmpfield=fieldnames(Ps);
        Ps.(tmpfield{3+namelist{jj}})=tmpval;
    else
        if(isempty(strfind(namelist{jj},'.')) && isempty(strfind(namelist{jj},'(')) && ~strcmp(namelist{jj},'Vs')) 
            Ps.(namelist{jj}) = tmpval;  % for parameters in Ps (Prefereable)
        else
            eval(sprintf('%s=tmpval;',namelist{jj}));   % Use eval, less safe/stable
        end;
	end;
end;


end
