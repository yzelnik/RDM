function bfs=StabBfByLP(bfs,pnts)
% Add fake stability values for a bifurcation array
% bfs=StabBfByLP(bfs,pnts)

if(nargin<2)
	pnts = 1;
end;

pnts=pnts(:)';
pnts = [1 sort(pnts) size(bfs,1)];

col = size(bfs,2)+1;
curstab=0;
for ii=1:length(pnts)-1
    bfs(pnts(ii):pnts(ii+1),col)=curstab;
    %disp([curstab pnts(ii:ii+1)])
    curstab=1-curstab;
end;

end
