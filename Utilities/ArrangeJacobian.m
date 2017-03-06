function Vs=ArrangeJacobian(Vs,Ps,Es)
% Arrange a Jacobian in different formats, according to Es.JacForm
% 0 = in a sparse matrix, with diagonals (same result if Es.JacForm doesn't exist)
% 1 = Just the data, with the first column for the derivative of all
%       equations by the first variable, looks like: [[UdU;VdU] [VdU; VdV]]


if(~isfield(Es,'JacForm') || Es.JacForm==0)
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
end;
    

end
