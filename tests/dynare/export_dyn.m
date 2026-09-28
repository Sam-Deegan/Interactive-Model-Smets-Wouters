load dyn_out.mat
f=fopen('dyn_irf.csv','w'); fprintf(f,'shock,variable,quarter,value\n');
v={'y','lab','pinf','r','c','inve','w','yf'}; s={'ea','eb','eg','eqs','em','epinf','ew'};
for i=1:7, for j=1:8, nm=[v{j} '_' s{i}];
 if isfield(oo_.irfs,nm), x=oo_.irfs.(nm); else x=zeros(1,20); end
 for t=1:20, fprintf(f,'%s,%s,%d,%.12e\n',s{i},v{j},t-1,x(t)); end
end, end
fclose(f);
c=oo_.conditional_variance_decomposition; size(c)
h=[1 2 3 4 5 6 10 11 40 41 100 101];
f=fopen('dyn_cvd.csv','w'); fprintf(f,'variable,horizon,shock,share\n');
exo=cellstr(M_.exo_names);
for j=1:8, for k=1:numel(h), for i=1:7
 fprintf(f,'%s,%d,%s,%.10f\n',v{j},h(k),exo{i},c(j,k,i)); end, end, end
fclose(f);
disp(exo')
