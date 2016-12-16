function [Vs,Ps,Es]=SaveParmList(Vs,Ps,Es,vallist,namelist)
% Save a list of parameters into Ps/Es
% [Vs,Ps,Es]=SavePrmList(Vs,Ps,Es,vallist,namelist)
% vallist is a vector of parameter values (otherwise Es.PrmList is used)
% namelist is a cell array of parameter names (otherwise Es.BfPrm is used)
% If Es.PrmInCell exists and is non-zero, then it points to values at Es.CellData
% A proper text name is assumed to be a Ps variable (unless it is "Vs")
% A numerical value acts as a pointer in the Ps structure (offset by +3)
% Anything else is directly evaluated using eval (not as safe) such as: Ps.Ds(2)
if(nargin<4) vallist = Es.PrmList; end;
if(nargin<5) namelist = Es.BfPrm; end;

% Do we have parameters who's values are given in cell-arrays?
Es=InsertDefaultValues(Es,'PrmInCell',zeros(length(namelist),1));

for jj=1:length(namelist)
    % Get value from vallist, one way or the other
	if(Es.PrmInCell(jj))
        tmpval=Es.CellData{Es.PrmInCell(jj)}{vallist(jj)};
    elseif(iscell(vallist))
        tmpval=vallist{jj};
    else
        tmpval=vallist(jj);
	end;
    % Update value into Ps/Es
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
