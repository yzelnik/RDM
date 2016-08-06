function bfs=MergePhases(bfs,cols,mults)
% merge more than one phase column into one
mults = [mults(:) ; ones(length(cols),1)];

for ii=1:size(bfs,1)
	for jj=1:size(bfs,2)
		temp = bfs{ii,jj};
		if(length(temp>0))
			endcol = temp(:,cols(1))*mults(1);
			for kk=2:length(cols)
				endcol = endcol + temp(:,cols(kk))*mults(kk);
			end;
			bfs{ii,jj} = [temp endcol];
		end;
	end;
end;
