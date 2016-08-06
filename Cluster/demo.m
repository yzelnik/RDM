function demo(par)

mkdir results;
cd results;

fid = fopen('demo_results.txt','a');
fprintf(fid,'%d\n',par);
fclose(fid);

end
