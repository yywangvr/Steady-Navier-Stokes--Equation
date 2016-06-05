function [C ,f]= ConvSource(X,T,referenceElement,velo)
% C = ConvectionMatrix(X,T,referenceElement,velo)
% Convection Matrix for a 2D Navier-Stokes problem
%
% X,T: nodal coordinates and connectivities for velocity
% referenceElement: reference element properties (quadrature, shape functions...)
% velo: velocity field


ngaus = referenceElement.ngaus;
wgp = referenceElement.GaussWeights;
N = referenceElement.N;
Nxi = referenceElement.Nxi;
Neta = referenceElement.Neta;

% reshape velocity field 
u_shape=[velo(1:2:end-1),velo(2:2:end)];

% Number of nodes and elements
nPt = size(X,1); 
nElem = size(T,1); 
nen=size(T,2);

C = zeros(nPt);
f = zeros(nPt,1);

% Loop on elements
for ielem = 1:nElem
    % Global number of the nodes in element ielem
    Te = T(ielem,:);
    % Coordinates of the nodes in element ielem
    Xe = X(Te,:);
    % Velocity at the element's nodes
    Conve = u_shape(Te,:);     
    % Element matrix
    [Ce,fe] = EleConvMatrix(Conve,Xe,nen,ngaus,wgp,N,Nxi,Neta);
    % Assemble the contrinbution of the element matrix
    C(Te, Te) = C(Te, Te) + Ce;
    f(Te)=f(Te)+fe;
end


end



function [Ce,fe] = EleConvMatrix(Conve,Xe,nen,ngaus,wgp,N,Nxi,Neta)
%

Ce = zeros(nen);
fe = zeros(nen,1);
% Loop on Gauss points 
for ig = 1:ngaus
   N_ig = N(ig,:);
    Nxi_ig = Nxi(ig,:);
    Neta_ig = Neta(ig,:);
    Jacob = [Nxi_ig*(Xe(:,1))	Nxi_ig*(Xe(:,2))
             Neta_ig*(Xe(:,1))	Neta_ig*(Xe(:,2))];
    dvolu = wgp(ig)*det(Jacob);
    res = Jacob\[Nxi_ig;Neta_ig];
    Nx = res(1,:);
    Ny = res(2,:);
    
     a = N_ig*Conve;
     ax = a(1); ay = a(2);
%     Pe_x=ax*h/(2*mu);
%     Pe_y=ax*h/(2*mu);
%     alpha_x = coth(Pe_x)-1/Pe_x;
%     alpha_y = coth(Pe_y)-1/Pe_y;
%     tau = h*(ax*alpha_x + ay*alpha_y)/2;
%     nubar_p = tau*norm(a)*norm(a);
    
    
    aGradN = (ax*Nx) + (ay*Ny);
    
    aux=N_ig*Conve;
    f_ig=update_Source(aux);
    
    Ce = Ce +( N_ig'*aGradN)*dvolu;
    fe = fe + N_ig'*(f_ig*dvolu);
  
end

end



function s=update_Source(aux)

%compute the norm in rows
nrm=norm(aux);

% compute updated source term
temp=1+exp(-10*(nrm-0.5));
s=1./temp;

end

