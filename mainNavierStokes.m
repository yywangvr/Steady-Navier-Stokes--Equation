% This program solves Stokes problem in a square domain


clear; close all; clc

addpath('Func_ReferenceElement')

dom = [0,2,0,3]; 


% Element type and interpolation degree
% (0: quadrilaterals, 1: triangles, 11: triangles with bubble function)
%elemV = 0; degreeV = 2; degreeP = 1;
 elemV = 1; degreeV = 2; degreeP = 1;
% elemV = 11; degreeV = 1;  degreeP = 1; 
if elemV == 11
    elemP = 1; 
else
    elemP = elemV; 
end
referenceElement = SetReferenceElementStokes(elemV,degreeV,elemP,degreeP); 
disp('');
hh = cinput('Spatial mesh size',0.2);
disp('');
nu=cinput('Diffusion coefficient nu for velocity',0.02);
disp('');
mu=cinput('Diffusion coefficient nu for velocity',1);

nx=(dom(2)-dom(1))/hh;
ny=(dom(4)-dom(3))/hh;
adapted = 0;
[X,T,XP,TP] = CreateMeshes(dom,nx,ny,referenceElement,adapted);

figure; PlotMesh(T,X,elemV,'b-');
figure; PlotMesh(TP,XP,elemP,'r-');
%% 

% Matrices arising from the discretization
[G,f,Kh] = FEM_system(X,T,XP,TP,referenceElement);
Kh=Kh*mu;
 
[ndofP,ndofV] = size(G); 

% Matrix and r.h.s vector to impose Dirichlet boundary conditions using
% Lagrange multipliers
[A_DirBC, b_DirBC, nDir, confined] = BC(X,dom,ndofV);
[Ah,bh,nDirh]=BC_h(X,dom);

if confined
   nunkP = ndofP-1;
   disp(' ')
   disp('Confined flow. Pressure on lower left corner is set to zero');
   G(1,:) = [];
else
   nunkP = ndofP;
end
nunkV = ndofV; 
btot = [f ; b_DirBC ; zeros(nunkP,1)];

% Initial solution
disp(' ')

lambda = zeros(nDir,1);
lambdah = zeros(nDirh,1);
pres = zeros(nunkP,1);
u = zeros(ndofV,1);
rho=ones(size(X,1),1);



iter = 0; tol = 0.5e-6; 

while iter < 100
    fprintf('Iteration = %d\n',iter);
    
    solh=[rho;lambdah];
    
    [Ch,fh]= ConvSource(X,T,referenceElement,u);
    
    Atoth=[Ch+Kh Ah';Ah zeros(nDirh)];
    ftoth = [fh;bh];

    res=ftoth-Atoth*solh;
     Atoth=sparse(Atoth);
     res=sparse(res);
     
     [Lh,Uh]=lu(Atoth);
    
    rhoInc=Uh\(Lh\res);
    
    rho=rho+rhoInc(1:length(fh));
    lambdah=lambdah+rhoInc(length(fh)+1:length(fh)+nDirh);
    
     K = Stiff_nu( X,T,referenceElement,rho,nu );
  
    Atot = [K        A_DirBC'             G'
            A_DirBC    zeros(nDir,nDir)     zeros(nDir,nunkP)
            G          zeros(nunkP,nDir)    zeros(nunkP,nunkP)];
    btot=[zeros(ndofV,1);b_DirBC; zeros(nunkP,1)];
    Atot=sparse(Atot);
    btot=sparse(btot);
    [L,U]=lu(Atot);
    
    sol=U\(L\btot);
    
    u=sol(1:ndofV);
    pres=sol(ndofV+nDir+1:ndofV+nDir+nunkP);
    
    
    
    delta1 = max(abs(rhoInc)); 
    delta2 = max(abs(res)); 
    
    if delta1 < tol && delta2 < tol
        fprintf('\nConvergence achieved in iteration number %g\n',iter); 
        break
    end
    
    iter=iter+1;
    
end

% Postprocess
if confined
    pres = [0; pres]; 
end

velo=[u(1:2:end-1),u(2:2:end)];


nPt = size(X,1); 
figure(3); 
quiver(X(1:nPt,1),X(1:nPt,2),velo(1:nPt,1),velo(1:nPt,2));
hold on 
plot(dom([1,2,2,1,1]),dom([3,3,4,4,3]),'k')
axis equal; axis tight

PlotStreamlines(X,velo,dom); 

pres=full(pres);
if degreeP == 0
    PlotResults(X,T,pres,referenceElement.elemP,referenceElement.degreeP)
else
    PlotResults(XP,TP,pres,referenceElement.elemP,referenceElement.degreeP)
end



xx  = reshape(X(:,1), degreeV*nx+1, degreeV*ny+1)'; 
yy  = reshape(X(:,2), degreeV*nx+1, degreeV*ny+1)';

%figure(3)


figure;
rhoo = reshape(rho, degreeV*nx+1, degreeV*ny+1)';  
surface(xx,yy,rhoo,'FaceColor','interp');
set(gca,'FontSize',16)
grid on
view(3)
pause(0.1)




