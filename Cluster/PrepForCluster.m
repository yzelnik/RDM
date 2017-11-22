function PrepForCluster(filename,varargin)
% Prepare a .mat file for loading during cluster runs
% First input is the filename to save, other inputs are for either matlab 
% files to run/load (e.g. Data.mat) or commands to run (e.g. ind = 3;)

% delete? Make sure the 3 basic variables exist, at the very least
Vs = [];
Ps = [];
Es = [];

varargin{nargin}=filename;  % adding filename to the end. Just a precaution

% Go over inputs
for ii=1:nargin-1
    input=varargin{ii};
    %disp(input)
    if(exist(input))  % This is a valid file name
        if(strcmpi(input(end-1:end),'.m')) % a matlab script
            eval(input(1:end-2)); % run the script
        elseif (strcmpi(input(end-3:end),'.mat')) % a matlab data file
            load(input);    % load the file
        else
            warning('Do not recognize file %s. Ignored it.',input);
        end;
    else
        %loc = strfind(input, '(');
        %if(~isempty(loc) && exist(input(1:loc-1)))
        eval([input ';']); % run this function (including parameters)
        %else
        %    loc = strfind(input, '=');
        %end;
    end;  
end;


filename=varargin{nargin};  % switching back (see begining of function)
if(~isfield(Es,'FileOut') || isempty(Es.FileOut))
    Es.FileOut = sprintf('clrun_%s',filename);  % Define output filename if needed
end;

% clear unnecessary variables from workspace, and save to file
clear('ii','input','varargin');
save(sprintf('%s_tmp.mat',filename));
%disp(sprintf('creating file: %s',filename));
end

