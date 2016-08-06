function [states, indx]=ReadAutoStates(filename,Ps,Es,pnts,varargin)
% Read AUTO states file (usually s.name or fort.8)
% points=ReadAutoStates(filename,Ps,Es,pnts)
% Returns a presentation of a set of states, identified by pnts

% Update online if necessary
[~,Ps,Es]=UpdateParameters([],Ps,Es,varargin{:});

if(~isfield(Es,'badtext')) % fix AUTO's bad-text problem (for very small/big values)
   Es.badtext=0;
end;

% Initialization
if(nargin<4) 
	pnts=0;
end;
abit    = 0.001;        % a very small value, used for interpolation
linelen = 16;           % the number of parms in an info line (before each state)
vrnum   = Ps.Vnum;      % how many variables to read per row
fsize   = Ps.Nx;        % how many points to interpolate into
maxpnt  = max(pnts);    % how many different states to try and read

% read file and read the info-line of the first state
fin   = fopen(filename,'r');
infoline  = fscanf(fin,'%f',linelen);

% get some general parameter from this infoline
rsz   = infoline(7);    % number of points (rows) saved per state
wid   = infoline(8);    % number of variables (columns)
totrowlen  = infoline(9)+1;  % how many rows per state in total - unused
extra = infoline(5)*2+infoline(12);    % number of extra numbers at the end 
lines = infoline(14)+infoline(15);     % number of lines this takes

tmpst = zeros(wid,rsz); % this will hold the data for each state we read

ind    = 0;
states = [];
indx   = [];
while((length(infoline)==linelen) && ((pnts(1)==0)||(ind<maxpnt)))
	indx = [indx infoline(2)]; % Add this state's index
    % read data of state
    if(Es.badtext==0)
    	[tmpst(:),~] = fscanf(fin,'%f',rsz*wid);    % Faster better method
    else
        tmpst = ReadSingleState(fin,rsz,wid)';  % Much slower, but works with bad files
    end;
	% read things for interpolation, and then do it
	Xs = [0-abit  tmpst(1,:)  1+abit ]';
	Ys = [tmpst(2:vrnum+1,1) tmpst(2:vrnum+1,:) tmpst(2:vrnum+1,end)]';
	states = cat(3,states,interp1(Xs,Ys,(0:(fsize-1))'/fsize));
	
    % read through data that comes at the end of each state, that we do not use
    if(Es.badtext==0)
    	 fscanf(fin,'%f',rsz*(wid-1)+extra);
    else
        for ii=1:(rsz+lines) fgetl(fin); end;
    end;
	
    % read next infoline of new state
    infoline  = fscanf(fin,'%f',linelen);

	ind = ind+1;
end;
% if not all points were asked for, then return only the specific list
if(pnts~=0)
	states = states(:,:,pnts);
end;

fclose(fin);

end

%%%%%%%%%%%%%%%%%% AUX FUNCTIONS %%%%%%%%%%%%%%%%%%

function st=ReadSingleState(fin,pointnum,fieldnum)
% This function attempts to read the state even in "difficult conditions"
numtxtsz=19;  % number of characters representing 1 number
offset=4;     % initial offset
tmpstr = fgetl(fin); % Clear the first line (half read from before)
for ii=1:pointnum
    tmpstr = fgetl(fin);
    for jj=1:fieldnum  % Go over fields in a single line
        substr = tmpstr((1:numtxtsz)+numtxtsz*(jj-1)+offset);
        if(sum(substr=='E')==0)
            substr(end-4)='E';  % Fix in the missing E in the text
        end;
        st(ii,jj)=sscanf(substr,'%f',1); % Make text into a (float) number
    end;
end;


end