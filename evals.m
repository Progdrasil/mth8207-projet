function [kQ,cQ,bQ,fQ]=evals(xq,k,c,b,f)
kQ=k(xq);
cQ=c(xq);
bQ=b(xq);
fQ=f(xq);