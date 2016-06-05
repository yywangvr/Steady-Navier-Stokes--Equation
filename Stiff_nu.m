function K = Stiff_nu( X,T,referenceElement,rho,nu )

elem = referenceElement.elemV;
ngaus = referenceElement.ngaus;
wgp = referenceElement.GaussWeights;
N = referenceElement.N;
Nxi = referenceElement.Nxi;
Neta = referenceElement.Neta;
ngeom = referenceElement.ngeom; 

[nElem,nenV] = size(T);
% Number of nodes
nPt_V = size(X,1);
if elem == 11
    nPt_V = nPt_V + nElem; 
end
nedofV = 2*nenV; 
ndofV = 2*nPt_V;
K = zeros(ndofV,ndofV);



% Loop on elements
for ielem = 1:nElem
    % Global number of the nodes in element ielem
    Te = T(ielem,:);
   
    % Coordinates of the nodes in element ielem
    Xe = X(Te(1:ngeom),:);
    rhoe=rho(Te(1:ngeom),:);
    % Degrees of freedom in element ielem
    Te_dof = reshape([2*Te-1; 2*Te],1,nedofV);
  
    
    % Element matrices
    [Ke] = EleMatStokes(Xe,ngeom,nedofV,ngaus,wgp,N,Nxi,Neta,rhoe,nu);
    
    % Assemble the element matrices
    K(Te_dof, Te_dof) = K(Te_dof, Te_dof) + Ke;
   
end


end

function [Ke] = EleMatStokes(Xe,ngeom,nedofV,ngaus,wgp,N,Nxi,Neta,rhoe,nu0)
%

Ke = zeros(nedofV,nedofV);

% Loop on Gauss points 
for ig = 1:ngaus
    N_ig    = N(ig,:);
    Nxi_ig  = Nxi(ig,:);
    Neta_ig = Neta(ig,:);
     
    Jacob = [
        Nxi_ig(1:ngeom)*(Xe(:,1))	Nxi_ig(1:ngeom)*(Xe(:,2))
        Neta_ig(1:ngeom)*(Xe(:,1))	Neta_ig(1:ngeom)*(Xe(:,2))
        ];
    dvolu = wgp(ig)*det(Jacob);
    res = Jacob\[Nxi_ig;Neta_ig];
    nx = res(1,:);
    ny = res(2,:);
    
    % Gradient
    Nx = [reshape([1;0]*nx,1,nedofV); reshape([0;1]*nx,1,nedofV)];
    Ny = [reshape([1;0]*ny,1,nedofV); reshape([0;1]*ny,1,nedofV)];
    
    aux=N_ig(1:ngeom)*rhoe;

    nu=update_nu(aux,nu0);
    Ke = Ke + nu*(Nx'*Nx+Ny'*Ny)*dvolu; 
 
end


end


function nu=update_nu(aux,nu0)

temp=1+exp(-10*(aux-0.5));
nu=nu0+nu0*(1./temp);

end