n = 500;
m = 1000;
A = complex(randn(n, n, m), randn(n, n, m));
tic
for i = 1:m
    A(:,:,i) * A(:,:,i);
end
toc
tic
pagemtimes(A, A);
toc