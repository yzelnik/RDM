function Es=InsertDefaultValues(Es,varargin)
% Add default values
% Es=InsertDefaultValues(Es,field1,value1,field2,value2,...)

% How many field+value sets do we have?
setnum=floor((nargin)/2);

% Run through all sets
for ii=1:setnum
     % Does field exist/is it empty?
    if(~isfield(Es,varargin{ii*2-1}) || isempty(Es.(varargin{ii*2-1})))
        % If so, update it with the default value
        Es.(varargin{ii*2-1})=varargin{ii*2}; 
    end;
end;


end