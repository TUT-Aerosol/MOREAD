% make a model coagulation sink for the model discrete

% v = dmpsload(datenum(1996,3,13));
v = load('sam20020414.sum');

V = Dlog_to_N(v);

dp = V(1,3:end);
N =  V(67,3:end).*0.5;

CO = Coags_dR(1e-9,dp,N,273)




fid = fopen('SINKDIST_05.TXT','w');
% return


fprintf(fid,'%i\n',length(dp))

for i = 1:length(dp)
	b1 = dp(i)./(10.^(floor(log10(dp(i)))));
	b2 = floor(log10(dp(i))); 
	if b2>=0,
		fprintf(fid,'%6.5fE+%i ',b1,b2)
	else
		fprintf(fid,'%6.5fE%i ',b1,b2)
	end
end
	fprintf(fid,'\b\n')
for i = 1:length(dp)
	b1 = N(i)./(10.^(floor(log10(N(i)))));
	b2 = floor(log10(N(i))); 
	if b2>=0,
		fprintf(fid,'%6.5fE+%i ',b1,b2)
	else
		fprintf(fid,'%6.5fE%i ',b1,b2)
	end
end
	fprintf(fid,'\b\n')

fclose(fid)
