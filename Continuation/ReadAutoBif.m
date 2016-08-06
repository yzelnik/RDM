function points=ReadAutoBif(filename,inds,mult,badtext)
% Read AUTO Bifurication file (usually b.name or fort.7)
% points=ReadAutoBif(filename,inds)
% Return the values of fields given in inds of all the points in the file
if(nargin<3) mult=0; end; % seperate to branches?
if(nargin<4) badtext=0; end; % try to repair bad-text?

% Initialization
if(nargin<2)
	inds=[1,2];
end;

% Open file
fin=fopen(filename);
tline=0;

% Build reading input format
inform = '%d %d %d %d ';
for ii=1:max(inds)
	inform = [inform '%f '];
end;

% Prepare for loop
tline=fgetl(fin);
points=[];
counter=[];
% Go over each line in file
while (tline~=-1)
	% Read which type of line this is
	flag = sscanf(tline,'%d/n');
    
	% If this line has info (first number not zero)
	if(flag~=0)
        if(badtext>0)  % Fix up the problems with the missing "E"
            tline=RepairLine(tline,max(inds)); 
        end;
		% Scan line and take in the fields specified in inds
		tmp = sscanf(tline,inform);
		points = [points;tmp(inds+4)'];
	end;
    counter = [counter;flag];
	% Read next line
	tline=fgetl(fin);
end;

if(mult)
    mult=max(counter);
    tmp=points;
    counter=nonzeros(counter); % only non-zeros values
    %size(counter)
    %size(tmp)
    points={};
    for ii=1:mult
        points{ii} = tmp((counter==ii),:);
    end;
    
end;

fclose(fin);

end

%%%%%%%%%%%%%%%%%% AUX FUNCTIONS %%%%%%%%%%%%%%%%%%

function newline=RepairLine(oldline,fieldnum)
% This function attempts to read the state even in "difficult conditions"
numtxtsz = 19;  % number of characters representing 1 number
offset   = 19;  % initial offset
place    = 4;   % where (from the end) should the E appear?  
newline  = oldline(1:offset);

for jj=1:fieldnum  % Go over fields in a single line
    substr = oldline((1:numtxtsz)+numtxtsz*(jj-1)+offset);
	if(sum(substr=='E')==0)
        substr(end-place)='E';  % Fix in the missing E in the text
    end;
    newline = [newline substr];
end;

end

