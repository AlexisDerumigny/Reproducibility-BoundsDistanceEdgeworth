
%df      = [5 5 5 5 5];
%lambda  = [1/5 1/5 1/5 1/5 1/5];

n = 5000;
df_original = 5;
df = df_original * ones(n, 1);
var = df_original / (df_original - 2);
lambda = (1/sqrt(var * n)) * ones(n, 1);
funtype = 1;
x=(-5):0.01:5;
[cdf,x] = tdist(x,df,lambda,funtype);

plot(x,cdf);

xstr = num2str(x,'%.3f');
cdfstr = num2str(cdf,'%.15f');
toWrite = [["x";xstr],["cdf";cdfstr]];

fileName = append("cdf_Student_n", num2str(n, '%.0f'), "_", ...
    num2str(df_original,'%.0f'), "df.csv");

writematrix(toWrite, fileName);


n = 5000;
df_original = 8;
df = df_original * ones(n, 1);
var = df_original / (df_original - 2);
lambda = (1/sqrt(var * n)) * ones(n, 1);
funtype = 1;
x=(-5):0.01:5;
[cdf,x] = tdist(x,df,lambda,funtype);

plot(x,cdf);

xstr = num2str(x,'%.3f');
cdfstr = num2str(cdf,'%.15f');
toWrite = [["x";xstr],["cdf";cdfstr]];

fileName = append("cdf_Student_n", num2str(n, '%.0f'), "_", ...
    num2str(df_original,'%.0f'), "df.csv");

writematrix(toWrite, fileName);
