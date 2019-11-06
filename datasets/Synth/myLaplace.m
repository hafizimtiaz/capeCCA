function sample = myLaplace(b,m,n)

% variance sigma^2 = 2*b^2
% m = number of rows in output sample
% n = number of cols in output sample
lambda=1/b;
siz=[m n];

zz = rand(siz);
xx = zeros(siz);
in = zz <=.5;
ip = zz > .5;
xx(in) =  1/lambda *log(2*zz(in));
xx(ip) = -1/lambda *log(2*(1-zz(ip)));

sample=xx;