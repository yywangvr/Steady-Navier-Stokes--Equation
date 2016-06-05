function [A, b, nDir, confined] = BC(X,dom,n)
% [A, b, nDir, confined] = BC(X,dom,n)
% Matrices to impose Dirichlet boundary conditions using Lagrange
% multipliers on a rectangular domain
% Input: 
%    X: nodal coordinates
%    dom: domain description [x1,x2,y1,y2]
%    n: number of velocity degrees of freedom
% Output:
%    A,b: matrix and r.h.s. vector to impose the boundary conditions using
%         Lagrange multipliers
%    nDir: number of prescribed degrees of freedom
%    confined: 

x1 = dom(1);  x2 = dom(2); 
y1 = dom(3);  y2 = dom(4);
%x2 = dom(2); y1 = dom(3);  y2 = dom(4); 

tol = 1e-6; 

nodesX1 = find(abs(X(:,1)-x1) < tol); 
nodesX2 = find(abs(X(:,1)-x2) < tol); 
nodesY1 = find(abs(X(:,2)-y1) < tol & abs(X(:,1)-x1) > tol & abs(X(:,1)-x2) > tol);
nodesY2 = find(abs(X(:,2)-y2) < tol & abs(X(:,1)-x1) > tol & abs(X(:,1)-x2) > tol); 



confined = 1; 
C = [
    2*nodesX1-1; 2*nodesX1
    2*nodesX2-1; 2*nodesX2
    2*nodesY1-1; 2*nodesY1
    2*nodesY2-1; 2*nodesY2
     ]; 
nDir = length(C);  
A = zeros(nDir,n); 
A(:,C) = eye(nDir); 
b = [
    zeros(size(nodesX1)); zeros(size(nodesX1))
    zeros(size(nodesX2)); -ones(size(nodesX2))
    zeros(size(nodesY1)); zeros(size(nodesY1))
    zeros(size(nodesY2));   zeros(size(nodesY2))
    ]; 

end




