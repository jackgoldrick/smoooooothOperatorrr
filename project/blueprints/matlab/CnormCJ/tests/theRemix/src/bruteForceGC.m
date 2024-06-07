%Effectively get norms plot as good as is possible, not caring about the
%computation cost
function [norms, vMax, p] bruteForceGC(cMatrix, pIn, dp, err)
    p = generatePValues(pIn);

    %Include code to create big set of p
    index2 = findIndex(p,2);
    
    pBig = p(index2:length(p));

    cMatrixP = conj(transpose(cMatrix));

    %Write for loop to solve through p to 2 for matrix and its transpose
    
    norms = zeros(1, length(p));
    vMax = zeros(length(cMatrix(1,:)),length(p(index2:end)));
    vMaxP = zeros(length(cMatrixP(1,:)),index2-1);
    
    [norms(index2-1), vMaxP(index2-1)] = genComparison(cMatrixP, p(index2-1)/(p(index2-1)-1), err, 50, 0);
    [norms(index2), vMax(1)] = genComparison(cMatrix, p(index2), err, 50, 0);

    for j = index2-2:-1:1
        q = p(j)/(p(j)-1);
        [norms(j), vMaxP(j)] = genComparison(cMatrix', q, err, 50, vMaxP(j+1));
    end
    for j = index2+1:length(p)
        [norms(j),vMax()]
end


%Upload generate p values to library

%Register for outlook on email
%Transfer drive files to outlook, or ensure they've been transferred
%Fix findIndex error


%%   for j = 1:index2-1
        %norms(j) = genComparison(cMatrix, p(j), err, 50, 0);
        %norms'(j) = genComparison(cMatrix', p(j), err, 50, 0);
    %end