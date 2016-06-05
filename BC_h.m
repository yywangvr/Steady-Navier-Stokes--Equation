function [ A,b,nDir ] = BC_h( X,dom)

x1 = dom(1);  x2 = dom(2); 
tol = 1e-6; 
 
%yM=1.5;
nodesX1 = find(abs(X(:,1)-x1) < tol ); 
nodesX2 = find(abs(X(:,1)-x2) < tol); 

C = [nodesX1, ones(size(nodesX1));
    nodesX2, zeros(size(nodesX2))];


nDir = size(C,1); 
nPt = size(X,1); 
A = zeros(nDir,nPt);
A(:,C(:,1)) = eye(nDir); 
b = C(:,2);



end

