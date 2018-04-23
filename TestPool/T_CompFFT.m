function comp = T_CompFFT(Vs,Ps,Es,varargin)
% Calculate the FFT-component(s) of a given state
% Where Es.CompFft contains the index(indices) of the component(s)
% Es.CompFft=0 is the uniform part, =1 is a single period, etc..

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if (~isfield(Es,'VarInd'))
    Es.VarInd = 1;
end;

% Preform the Fourier Transform
tmpfft=fft(Vs(:,Es.VarInd(1),1));

% Go over each component requested
for ii=1:length(Es.CompFft)
    if(Es.CompFft(ii)==0)   % for the zero-order component
        comp(ii)=mean(Vs(:,Es.VarInd(1),1));
    else                    % or any other
        comp(ii)=sum(abs(tmpfft([Es.CompFft(ii)+1 length(tmpfft)-Es.CompFft(ii)+1])))/length(tmpfft);
    end;
    
end;

end