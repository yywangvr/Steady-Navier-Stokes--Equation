function Element = SetReferenceElement(elem,p)
% Element = SetReferenceElement(elem,p)
% INPUT: 
%    elem: type of element (0: quadrilatera, 1: triangles)
%    p: interpolation degree
% OUTPUT: 
%    Element: struture with the reference element's properties

if elem == 0
    if p == 1
        ngaus = 4;
        Xe_ref = [-1,-1; 1,-1; 1,1; -1,1];
    elseif p == 2
        ngaus = 9; 
        Xe_ref =  [-1,-1; 1,-1; 1,1; -1,1; ...
            0,-1; 1,0; 0,1; -1,0; 0,0];
    else
        error('not available interpolation degree');
    end
elseif elem == 1
    if p == 1
        ngaus = 3; 
        Xe_ref = [1,0; 0,1; 0,0];
    elseif p == 2
        ngaus = 6; 
        Xe_ref = [1,0; 0,1; 0,0; 0.5,0.5; 0,0.5; 0.5,0];
    else
        error('not available interpolation degree');
    end
elseif elem == 11
    if p == 1
        ngaus = 3; 
        Xe_ref = [0,0; 1,0; 0,1];
    else
        error('not available interpolation degree');
    end
else
    error('unavailable element')
end

% Gauss points and weights
[zgp,wgp] = Quadrature(elem,ngaus);

% Shape functions and their derivatives (in reference coordinates)
[N,Nxi,Neta] = ShapeFunc(elem,p,zgp);

Element.degree = p;
Element.elem = elem;
Element.Xe_ref = Xe_ref;
Element.nen = size(Xe_ref,1); 
Element.ngaus = ngaus; 
Element.GaussPoints = zgp; 
Element.GaussWeights = wgp; 
Element.N = N; 
Element.Nxi = Nxi; 
Element.Neta = Neta; 

