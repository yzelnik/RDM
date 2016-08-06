function tscale = C_CalcConvergenceTimescale(Input,Ps,Es,varargin)
% Estimate the time scale of convergence (assuming exponential decay)

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;

if(~isfield(Es,'BfFields'))
 	Es.BfFields=[1,2];
end;

if(size(Input,3)>1)     % If we get state data, run the test on each state to form a history
    bfhist(1:size(Input,3)) = 1:size(Input,3);
    for ii=1:size(Input,3)
        temp = Es.TestFunc(Input(:,:,ii),Ps,Es);
        bfhist(ii,2:1+length(temp))=temp;
    end;
else                    % Or, assume we got a history
    bfhist = Input;
end;

% Find the max difference over time, assume convergence starts after it
[~,mind] = max(bfhist(:,Es.BfFields(2)));		
mind = max([mind round(size(bfhist,1)/2)]);     % add-hoc correction, don't use first half of data

% fit to a linear function (with the log of the values in the y axis)
coeffs   = polyfit(bfhist(mind:end,Es.BfFields(1)),log(bfhist(mind:end,Es.BfFields(2))),1);
% use the slope (of a logarithmic function to begin with) to get the timescale and also put in the overall time
tscale   = [-1/coeffs(1) bfhist(end,Es.BfFields(1))];		

%plot(bfhist(mind:end,1),log(bfhist(mind:end,fieldind)))

end