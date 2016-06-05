function [N,Nxi,Neta] = ShapeFunc(elem,p,z) 
% [N,Nxi,Neta] = ShapeFunc(elem,z)  
% N, Nxi, Neta: matrices storing the values of the shape functions on the Gauss points
%               of the reference element
%               Each row concerns to a Gauss point
% z: coordinates of Gauss points in the reference element
% elem:     type of element (0: quadrilatera, 1: triangles, 11: triangle+bubble function)
% p: interpolation degree

xi = z(:,1); eta = z(:,2); 

if elem == 0 % quadrilatera
    if p == 0
        N    = ones(size(xi));
        Nxi  = zeros(size(xi));
        Neta = zeros(size(xi));        
    elseif p == 1
        N    = [(1-xi).*(1-eta)/4, (1+xi).*(1-eta)/4, (1+xi).*(1+eta)/4, (1-xi).*(1+eta)/4]; 
        Nxi  = [(eta-1)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4]; 
        Neta = [(xi-1)/4, -(1+xi)/4,   (1+xi)/4,  (1-xi)/4 ]; 
    elseif p == 2
        N    = [xi.*(xi-1).*eta.*(eta-1)/4, xi.*(xi+1).*eta.*(eta-1)/4, ...
            xi.*(xi+1).*eta.*(eta+1)/4, xi.*(xi-1).*eta.*(eta+1)/4, ...
            (1-xi.^2).*eta.*(eta-1)/2,  xi.*(xi+1).*(1-eta.^2)/2,   ...
            (1-xi.^2).*eta.*(eta+1)/2,  xi.*(xi-1).*(1-eta.^2)/2,   ...
            (1-xi.^2).*(1-eta.^2)];
        Nxi  = [(xi-1/2).*eta.*(eta-1)/2,   (xi+1/2).*eta.*(eta-1)/2, ...
            (xi+1/2).*eta.*(eta+1)/2,   (xi-1/2).*eta.*(eta+1)/2, ...
            -xi.*eta.*(eta-1),          (xi+1/2).*(1-eta.^2),   ...
            -xi.*eta.*(eta+1),          (xi-1/2).*(1-eta.^2),   ...
            -2*xi.*(1-eta.^2)];
        Neta = [xi.*(xi-1).*(eta-1/2)/2,    xi.*(xi+1).*(eta-1/2)/2, ...
            xi.*(xi+1).*(eta+1/2)/2,    xi.*(xi-1).*(eta+1/2)/2, ...
            (1-xi.^2).*(eta-1/2),       xi.*(xi+1).*(-eta),   ...
            (1-xi.^2).*(eta+1/2),       xi.*(xi-1).*(-eta),   ...
            (1-xi.^2).*(-2*eta)];        
    else
        error('not available interpolation degree')
    end
elseif elem == 1 % triangle
    if p == 0
        N    = ones(size(xi));
        Nxi  = zeros(size(xi));
        Neta = zeros(size(xi));
    elseif p == 1
        N    = [xi,              eta,             1-xi-eta]; 
        Nxi  = [ones(size(xi)),  zeros(size(xi)), -ones(size(xi))];
        Neta = [zeros(size(xi)), ones(size(xi)),  -ones(size(xi))];
    elseif p == 2
        N    = [xi.*(2*xi-1),eta.*(2*eta-1),(1-2*(xi+eta)).*(1-(xi+eta)),4*xi.*eta,4*eta.*(1-(xi+eta)),4*xi.*(1-(xi+eta))]; 
        Nxi  = [4*xi-1,zeros(size(xi)),-3+4*(xi+eta),4*eta,-4*eta,4*(1-2*xi-eta)]; 
        Neta = [zeros(size(xi)),4*eta-1,-3+4*(xi+eta),4*xi,4*(1-xi-2*eta),-4*xi];         
    else
        error('not available interpolation degree')
    end   
elseif elem == 11 % triangle + bubble function
    if p == 1
        N    = [xi,eta,1-(xi+eta),27*xi.*eta.*(1-xi-eta)]; 
        Nxi  = [ones(size(xi)),zeros(size(xi)),-ones(size(xi)),27*eta.*(1-2*xi-eta)]; 
        Neta = [zeros(size(xi)),ones(size(xi)),-ones(size(xi)),27*xi.*(1-2*eta-xi)]; 
    else
        error('not available interpolation degree');
    end
    
else
    error('not available element')
end
