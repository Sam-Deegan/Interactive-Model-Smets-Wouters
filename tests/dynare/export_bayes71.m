load dyn_est71.mat
f=fopen('bayes_irf71.csv','w'); fprintf(f,'shock,variable,quarter,mean,median\n');
v={'y','lab','pinf','r'}; s={'ea','eb','eg','eqs','em','epinf','ew'};
for i=1:7, for j=1:4, nm=[v{j} '_' s{i}];
 if isfield(oo_.PosteriorIRF.dsge.Mean,nm), x=oo_.PosteriorIRF.dsge.Mean.(nm); m=oo_.PosteriorIRF.dsge.Median.(nm); else x=zeros(1,20); m=x; end
 for t=1:20, fprintf(f,'%s,%s,%d,%.8e,%.8e\n',s{i},v{j},t-1,x(t),m(t)); end
end, end
fclose(f);
disp(oo_.posterior_mean.shocks_std)
