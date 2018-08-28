function MakeStVideo(Vs,Ps,Es,varargin)
% Build a video animation from state-data (Vs) and write to file: Es.FileOut
% MakeStMovie(Vs,Ps,Es)
% Use Es.TitlesFrames and Es.TitlesText for titles along video
% Use Es.VideoSpeed and Es.VideoQuality to augment video

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

Es=InsertDefaultValues(Es,'PlotFunc',@plotst);

if(~isfield(Es,'VideoSpeed') || Es.VideoSpeed==0) % setup video speed?
	Es.VideoSpeed=1;
end;

FileName = Es.FileOut; % file name to write video to
FrameRate = Es.VideoSpeed; % define frame rate

% Setup title timing if necessary
if(isfield(Es,'TitlesFrames'))
    if(length(Es.TitlesFrames)<size(Vs,3)) % if only timing of title-changing is given
        for ii=1:length(Es.TitlesFrames)-1 % go through each title
            temp(Es.TitlesFrames(ii):Es.TitlesFrames(ii+1)-1)=ii;
        end;
        temp(Es.TitlesFrames(ii+1):size(Vs,3))=ii+1; % last title lasts to the end
        Es.TitlesFrames=temp; 
    end;
else
    Es.TitlesFrames=zeros(size(Vs,3),1); % no titles
end;

if(~sum(FileName=='.')) % If not file type was specified, use avi
    FileName = [FileName '.avi'];
end;

folder = pwd; % deine the folder
if(prod(FileName(end-2:end)=='mp4'))  % setup video writer
    writerObj = VideoWriter([folder '/' FileName],'MPEG-4');
else
    writerObj = VideoWriter([folder '/' FileName]);
end;
writerObj.FrameRate = FrameRate;
if(isfield(Es,'VideoQuality') && ~(Es.VideoQuality==0))
	writerObj.Quality = Es.VideoQuality; % define video quality if needed
end;

% start the writing and get figure
open(writerObj);
fig=gcf;

% Go over each frame
for k=1:size(Vs,3)
    % plot the state
    Es.PlotFunc(Vs(:,:,k),Ps,Es);
    if(Es.TitlesFrames(k)) % plot title if relevant
        title(Es.TitlesText{Es.TitlesFrames(k)});
    end;
    % get frame and write it
    frame = getframe(fig);
    writeVideo(writerObj,frame);
end

% finish-up writing
close(writerObj);

end
