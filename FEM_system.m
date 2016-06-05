function [G,f,Kh] = FEM_system(X,T,XP,TP,referenceElement)
% [K,G,f] = Stokes_system(X,T,XP,TP,referenceElement)
% Matrices K, G and r.h.s vector f obtained after discretizing a Stokes problem
%
% X,T: nodal coordinates and connectivities for velocity
% XP,TP: nodal coordinates and connectivities for pressure
% referenceElement: reference element properties (quadrature, shape functions...)


elem = referenceElement.elemV;
ngaus = referenceElement.ngaus;
wgp = referenceElement.GaussWeights;
N = referenceElement.N;
Nxi = referenceElement.Nxi;
Neta = referenceElement.Neta;
NP = referenceElement.NP; 
ngeom = referenceElement.ngeom; 

% Number of elements and number of nodes in each element
[nElem,nenV] = size(T);
nenP = size(TP,2); 

% Number of nodes
nPt_V = size(X,1);
if elem == 11
    nPt_V = nPt_V + nElem; 
end
nPt_P = size(XP,1);

% Number of degrees of freedom 
nedofV = 2*nenV; 
nedofP = nenP;
ndofV = 2*nPt_V; 
ndofP = nPt_P; 

%K = zeros(ndofV,ndofV);
G = zeros(ndofP,ndofV); 
f = zeros(ndofV,1);
Kh= zeros(nPt_V,nPt_V);
% Loop on elements
for ielem = 1:nElem
    % Global number of the nodes in element ielem
    Te = T(ielem,:);
    TPe = TP(ielem,:); 
    % Coordinates of the nodes in element ielem
    Xe = X(Te(1:ngeom),:);
    % Degrees of freedom in element ielem
    Te_dof = reshape([2*Te-1; 2*Te],1,nedofV);
    TPe_dof = TPe; 
    
    % Element matrices
    [Ge,Khe,fe] = EleMatStokes(Xe,ngeom,nenV,nedofV,nedofP,ngaus,wgp,N,Nxi,Neta,NP);
    
    % Assemble the element matrices
    %K(Te_dof, Te_dof) = K(Te_dof, Te_dof) + Ke;
    G(TPe_dof,Te_dof) = G(TPe_dof,Te_dof) + Ge; 
    f(Te_dof) = f(Te_dof) + fe;
    Kh(Te,Te)=Kh(Te,Te)+Khe;
end






function [Ge,Khe,fe] = EleMatStokes(Xe,ngeom,nenV,nedofV,nedofP,ngaus,wgp,N,Nxi,Neta,NP)
%

%Ke = zeros(nedofV,nedofV);
Ge = zeros(nedofP,nedofV);
fe = zeros(nedofV,1);
Khe= zeros(nenV,nenV);

% Loop on Gauss points 
for ig = 1:ngaus
    N_ig    = N(ig,:);
    Nxi_ig  = Nxi(ig,:);
    Neta_ig = Neta(ig,:);
    NP_ig = NP(ig,:); 
    Jacob = [
        Nxi_ig(1:ngeom)*(Xe(:,1))	Nxi_ig(1:ngeom)*(Xe(:,2))
        Neta_ig(1:ngeom)*(Xe(:,1))	Neta_ig(1:ngeom)*(Xe(:,2))
        ];
    dvolu = wgp(ig)*det(Jacob);
    res = Jacob\[Nxi_ig;Neta_ig];
    nx = res(1,:);
    ny = res(2,:);
    
	Ngp = [reshape([1;0]*N_ig,1,nedofV); reshape([0;1]*N_ig,1,nedofV)];
    % Gradient
   % Nx = [reshape([1;0]*nx,1,nedofV); reshape([0;1]*nx,1,nedofV)];
   % Ny = [reshape([1;0]*ny,1,nedofV); reshape([0;1]*ny,1,nedofV)];
    % Divergence
    dN = reshape(res,1,nedofV);

  
    Khe = Khe+(nx'*nx+ny'*ny)*dvolu;
    Ge = Ge - NP_ig'*dN*dvolu; 
    x_ig = N_ig(1:ngeom)*Xe; 
    f_igaus = SourceTerm(x_ig); 
    fe = fe + Ngp'*f_igaus*dvolu; 
end

