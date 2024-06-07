

A = randn(4096);
tic
w = inv(A);
toc

i = eye(4096);
tic
v = A \ i;
toc

l = max(abs(w-v), [], 'all')
