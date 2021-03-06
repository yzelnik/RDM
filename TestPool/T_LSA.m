function [res,evs,iternum]=T_LSA(Vs,Ps,Es,varargin)
% Linear Stability Analysis of a numerical solution
% [res,evs,ind]=T_LSA(Vs,Ps,Es)
% Uses the Ps.LocFunc and Ps.SpaFunc to construct a Jacobian
% Returns 1 if the state is stable, 0 otherwise (-1 if unsuccessful)
% evs returns first Es.EigNum eigenvalues (default #=4)
% iternum returns the number of iterations required for the analysis
% Use Es.AvoidErrors>0 to run the eigs function even if errors occur
% A value between 0 and 1 implies the max error (tolerance) allowed
% Otherwise, the value is the max number of repetitions to be attempted

% Update online if necessary
if(nargin>3) [Vs,Ps,Es]=UpdateParameters(Vs,Ps,Es,varargin{:}); end;


Es.LsaThresh=[Es.LsaThresh(:);0]; % pad with zero at the end
tmp=getjac(Vs,Ps,Es);
        %spy(tmp)
ind = 0;         % Used for tracking the number of calculating attempts
eignum = 4;
if( (isfield(Es,'EigNum')) && Es.EigNum)
    eignum=Es.EigNum;
end;

if(~isfield(Es,'JacNum'))
    Es.JacNum=0;
end;

% If Matlab crashes (because of eigs) - we can try to deal with this:
if( (isfield(Es,'AvoidErrors')) && Es.AvoidErrors)
	flag = 1;
	tol = eps;
	% If Es.AvoidErrors is an integer, repeat up to Es.AvoidErrors times
	if(~mod(Es.AvoidErrors,1))	
		itermax = Es.AvoidErrors;
		tolchange = 1;
	else	% If This is not an integer, than we want to change opts.tol with iterations, up to Es.AvoidErrors
		tolchange = 10;
		%tolmax = Es.AvoidErrors;
		itermax = ceil(log10(Es.AvoidErrors/tol));
	end;
	% Iterate running eigs, catching the error if it occurs
	while flag && ind<itermax
		ind=ind+1;
		flag=0;
		opts.tol = tol;
        
		try
            evs = eigs(CalculateJacobian(Vs,Ps,Es),eignum,1,opts);
           % if(~Es.JacNum)  % Use analytical jacobian?
    			%evs=eigs(Ps.LocFunc(Vs,Ps,Es)+Ps.SpaFunc(Vs,Ps,Es),eignum,1,opts);
            %    evs=eigs(getjac(Vs,Ps,Es),eignum,1,opts);
            %else            % or numeric one?
            %    evs=eigs(NumericJacobian(Vs,Ps,Es),eignum,1,opts);
            %end;
        catch
            flag = 1;
			evs = [0 0 0 0];
		end;
		tol = tol * tolchange;
	end;
	if(ind>1)	% If more than one iteration was required, give some output, either to the screen or file
			tmp = sprintf('*** Done with %d runs! Nx: %.4f, Lx: %.4f, tolerance: %.4e\n',ind,Ps.Nx,Ps.Lx,tol);
			if( (isfield(Es,'BfOut')) & Es.BfOut)
				fin=fopen(Es.BfOut,'a');	% Print to file named Es.BfOut
				fprintf(fin,tmp);
				fclose(fin);
			else	% Just print to screen
				fprintf(tmp);
			end;
	end;
    
elseif( (isfield(Es,'UseEig')) && Es.UseEig)
    % use eig instead of eigs
    tmp = eig(full(CalculateJacobian(Vs,Ps,Es)));
    tmp = sort(real(tmp),'descend');
    evs = tmp(1:eignum);
else
    % Notmal run of eigs
    evs = eigs(CalculateJacobian(Vs,Ps,Es),eignum,1);
%	if(~Es.JacNum)  % Use analytical jacobian?
 %       tmp=CalculateJacobian(Vs,Ps,Es);
        %spy(tmp)
  %      evs=eigs(tmp,eignum,1);
        %evs=eigs(Ps.LocFunc(Vs,Ps,Es)+Ps.SpaFunc(Vs,Ps,Es),eignum,1);
   %     disp(999)
%	else            % or numeric one?
%        evs=eigs(NumericJacobian(Vs,Ps,Es),eignum,1);%
%	end;
end;


if(sum(abs(evs))==0)	% If the run was not successful, than return -1
	res = -1;
else	% eigs was successful, so return either 1 (stable) or 0 (unstable), using the second eigenvalue
	temp = sort(real(evs));
    [~,zeromode] = min(abs(temp));
    if(temp(zeromode)<abs(Es.LsaThresh(1)))
        %disp(sprintf('found mode at %d, with %f',zeromode,temp(zeromode)));
        temp(zeromode)=-inf;    % Ignore zero mode ( zero eigvalue)
    end;
    
    res = (1-(max(temp)>Es.LsaThresh(2)))>0;
	%res = temp(end-1)<Es.LsaThresh;
end;

iternum = ind;
 

end

