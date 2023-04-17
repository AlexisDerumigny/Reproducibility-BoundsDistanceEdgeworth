n = 40000;
df_original = 5;

df = df_original * ones(n, 1);
var = df_original / (df_original - 2);

lambda = (1/sqrt(var * n)) * ones(n, 1);
funtype = 1;
x=(-5):0.1:5;

[cdf,x] = tdist(x,df,lambda,funtype);

plot(x,cdf);

toWrite = [["x";x],["cdf";cdf]];
writematrix(toWrite, "cdf_Student.csv");
