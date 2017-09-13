function Vs=ArrangeJacobian(Vs,Ps,Es)
% Switch between format of Jacobian:
% 1) in a sparse matrix, with diagonals 
% 2) just the data, with the first column for the derivative of all
%    equations by the first variable, looks like: [[UdU;VdU] [VdU; VdV]]

% determine which format the jacobian is currently in
if(size(Vs,1)>size(Vs,2))
    sparseform=0; % should be in form 2
elseif(size(Vs,1)==size(Vs,2))
    if(Ps.VarNum*Ps.Nx*Ps.Ny==size(Vs,1))
        sparseform=1; % it is in form 1
    else
        sparseform=0; % should be in form 2
    end;
else
    error('Jacobian appears to have the wrong format'); % nothing makes sense
end;


if(~sparseform) % switch from form 2 to form 1
	% write in a large sparse matrix format (along diagonals)
    num   = size(Vs,2);
    len   = size(Vs,1)/num;
    temp  = zeros(len*num,num*2-1);
    
    % Go over diagonals, and for each one go over each equation part
    for ii=1-num:num-1
        for jj=1:num
            %[jj jj-ii ii+num jj len] 
            if((jj-ii)>0 && (jj-ii)<=num)
                temp((1:len)+(jj-1)*len,ii+num)=Vs((1:len)+(jj-ii-1)*len,jj);
            end;
        end;
    end; 
    % Create the sparse matrix using the diagonals 
    Vs = spdiags(temp,(1-num:num-1)*len,len*num,len*num);
else  % switch from form 1 to form 2
    len=size(Vs,1)/Ps.VarNum;
    %tmp=spdiags(Vs);
    tmp=Vs;
    %imagesc(tmp>0)
    Vs=zeros(Ps.VarNum*len,Ps.VarNum);
    for ii=1:Ps.VarNum
        for jj=1:Ps.VarNum
            
            Vs((1:len)+(ii-1)*len,jj) = diag(tmp((1:len)+(ii-1)*len,(1:len)+(jj-1)*len));
            %[ii jj len]
            %Vs((1:len)+(ii-1)*len,jj) = tmp((1:len)+(jj-1)*len,jj+Ps.VarNum-ii);
        end;
    end;
end;
    
end
