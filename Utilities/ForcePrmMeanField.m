function Ps=ForcePrmMeanField(Ps)
% Force a Ps structure to have mean-field parameters (if any are not)

% Assume parameters start from the 4rd field in the Ps structure
defstart=4;
% Read parameters into an array
fn=fieldnames(Ps);
ind = defstart;
len=Ps.Nx*Ps.Ny;

% go through all parameters, until end of "hitting" Ps.Ds
while (ind<length(fn)) && ~strcmp(fn{ind},'Ds')
    % do we have a non-scalar parameter?
    if(size(Ps.(fn{ind}),1)==len) 
        Ps.(fn{ind}) = mean(Ps.(fn{ind}),1); % average it out to mean-field
        disp(Ps.(fn{ind}));
    end;
    ind=ind+1;
end;

end
