function [X,T,XP,TP] = CreateMeshes(dom,nx,ny,referenceElement,adapted)
% Uniform meshes in a rectangular domain
% Input:    
%   dom = [x1,x2,y1,y2]:  vertices' coordinates
%   nx,ny: number of elements in each direction
%   referenceElement: reference element's properties
%   adapted = 1 for a mesh which is refined near the boundary
% Output:   
%   X,T:   nodal coordinates and connectivities of the velocity mesh
%   XP,TP: nodal coordinates and connectivities of the pressure mesh

elemV = referenceElement.elemV; 
degreeV = referenceElement.degreeV; 
elemP = referenceElement.elemP;
degreeP = referenceElement.degreeP; 

if adapted == 1
    [X,T] = CreateAdaptedMesh(dom,nx,ny,elemV,degreeV); 
else
    [X,T] = CreateUniformMesh(dom,nx,ny,elemV,degreeV); 
end

if degreeP == 0
    nElem = size(T,1); 
    TP = (1:nElem)'; 
    XP = zeros(nElem,2); 
    for i = 1:nElem
        Te = T(i,:); 
        Xe = X(Te,:); 
        XP(i,:) = [mean(Xe(:,1)), mean(Xe(:,2))];
    end
elseif elemV == 11
    warning('only linear elements')
    XP = X; 
    TP = T(:,1:3); 
else
    if adapted == 1
        [XP,TP] = CreateAdaptedMesh(dom,nx,ny,elemP,degreeP); 
    else
        [XP,TP] = CreateUniformMesh(dom,nx,ny,elemP,degreeP); 
    end
end







