function PrintDimensionalization(Ps,NewParDef,VarNames)
% Print out the dimensionalization info from a matrix
if(nargin<3)
    for ii=1:Ps.VarNum
        VarNames{ii}=sprintf('v%d',ii);
    end;
end;
% Assume parameters start from the 3rd field in the Ps structure
defstart=3;
plen = size(NewParDef,2);

% Read parameters into an array
tempparms=[];
fn  = fieldnames(Ps);
ind = 1;         % follows parameter value (relevant if parameter holds more than one)
ii  = defstart;  % follows parameter name
while ind<plen
    parmsize = length(Ps.(fn{ii})(:));
    if(parmsize==1)
        parmnames{ind} = fn{ii};
        ind=ind+1;
    else
        for jj=1:parmsize
            parmnames{ind}= sprintf('%s_%d',fn{ii},jj);
            ind=ind+1;
        end;
    end;
    ii=ii+1;
end;


for ii=1:plen
    if(~isempty(nonzeros(NewParDef(ii,:))))
        fprintf('%s_new = \t',parmnames{ii});
        PrintRelations(parmnames,NewParDef(ii,:))
     end;
end;

for ii=1:Ps.VarNum
    fprintf('%s_new = \t',VarNames{ii});
    PrintRelations([VarNames(ii),parmnames],[1 NewParDef(plen+ii,:)])
    %fprintf('%s_new = %s * ',VarNames{ii},VarNames{ii});
    %PrintRelations(parmnames,NewParDef(plen+ii,:))
end;


   % PrintRelations([VarNames(ii),parmnames],[1 NewParDef(plen+ii,:)])
fprintf('x_new = \t');
PrintRelations([{'x'},parmnames],[1 NewParDef(plen+Ps.VarNum+1,:)])
fprintf('t_new = \t');
PrintRelations([{'t'},parmnames],[1 NewParDef(plen+Ps.VarNum+2,:)])
%(parmnames,NewParDef(plen+Ps.VarNum+2,:))

end


function PrintRelations(parmnames,relations)
    linearplus = find(relations==1);
    linearminus = find(relations==-1);
    
    if(~isempty(linearplus))    % Positive linear relations (par^1)
        fprintf('%s',parmnames{linearplus(1)});
        for ii=2:length(linearplus) 
            %disp([1 ii 1])
            fprintf('*%s',parmnames{linearplus(ii)});
        end;
    else
        if(~isempty(linearminus))
            fprintf('1');
        end;
    end;
    if(~isempty(linearminus))   % Negative linear relations (par^-1)
        fprintf(' / %s',parmnames{linearminus(1)});
        for ii=2:length(linearminus)
            %disp([0 ii 0])
            fprintf('*%s',parmnames{linearminus(ii)});
        end;
    end;
    % Take care of non-linear relations
    linearlen = length(linearplus) + length(linearminus);
    relations(linearplus)=0;
    relations(linearminus)=0;
    leftovers=find(relations);
    if(~isempty(leftovers))
        if(linearlen==0)
            fprintf('*--');
        end;
        for ii=1:length(leftovers)
            fprintf(' * %s^%.1f',parmnames{leftovers(ii)},relations(leftovers(ii)));
        end;
        
    end;
    
    fprintf('\n');
end

