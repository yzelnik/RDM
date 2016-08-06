function vallist=LoadParmList(Vs,Ps,Es,namelist)
% Loads a list of parameters from Ps/Es into vallist
% vallist=LoadParmList(Vs,Ps,Es,namelist)
% namelist is a cell array of parameter names (otherwise Es.BFpar is used)
% A proper text name is assumed to be a Ps variable (unless it is "Vs")
% A numerical value acts as a pointer in the Ps structure (offset by +3)
% Anything else is directly evaluated using eval (not as safe) such as: Ds(2)
if(nargin<4) namelist = Es.BFpar; end;

for jj=1:length(namelist)
	if(isnumeric(namelist{jj})) % Allow access to model-parameters by index
        tmpfield = fieldnames(Ps);
        tmpval = Ps.(tmpfield{3+namelist{jj}});
    else
        if(isempty(strfind(namelist{jj},'.')) && isempty(strfind(namelist{jj},'(')) && ~strcmp(namelist{jj},'Vs')) 
            tmpval = Ps.(namelist{jj});  % for parameters in Ps (Prefereable)
        else
            eval(sprintf('tmpval=%s;',namelist{jj}));   % Use eval, less safe/stable
        end;
	end;
    vallist(jj)=tmpval;
end;

end