function PlotResults(X,T,sol,elem,degree)

if elem == 1 && degree == 1
    figure; 
    trisurf(T,X(:,1),X(:,2),sol,'FaceColor','interp');
elseif degree == 0
    nElem = size(T,1); 
    if elem == 0
        nen = 4; 
    else
        nen = 3; 
    end
    figure; hold on
    for ielem = 1:nElem
        Te = T(ielem,1:nen); 
        Xe = X(Te,:); 
        zz = sol(ielem)*ones(nen,1);
        patch(Xe(:,1),Xe(:,2),zz,zz)
    end
else
    [nElem,nen] = size(T);
    if elem == 11
        ngeom = nen-1; 
    else
        ngeom = nen; 
    end
    if elem == 0
        npt = 2*degree+1; 
        x = linspace(-1,1,npt); 
        [x,y] = meshgrid(x,x); 
        pts = [reshape(x,npt^2,1), reshape(y,npt^2,1)]; 
        ptsEdge = [
            linspace(-1,1,npt)'            -ones(npt,1)
            ones(npt,1)                    linspace(-1,1,npt)'       
            flipud(linspace(-1,1,npt)')    ones(npt,1)
            -ones(npt,1)                   flipud(linspace(-1,1,npt)')
            ]; 
    else
        npt = 2*degree+1; 
        x = linspace(0,1,npt); 
        [x,y] = meshgrid(x,x); 
        pts = [reshape(x,npt^2,1), reshape(y,npt^2,1)]; 
        ind = pts(:,2) <= 1 - pts(:,1); 
        pts = pts(ind,:); 
        ptsEdge = [
            linspace(0,1,npt)'            1-linspace(0,1,npt)'
            flipud(linspace(0,1,npt)')    zeros(npt,1)
            zeros(npt,1)                  linspace(0,1,npt)'
            ];
    end
    tri = delaunay(pts(:,1), pts(:,2)); 
    N = ShapeFunc(elem,degree,pts);
    NEdge = ShapeFunc(elem,degree,ptsEdge);
    figure; hold on; 
    for ielem = 1:nElem
        Te = T(ielem,:); 
        Xe = X(Te(1:ngeom),:); 
        sol_e = sol(Te);
        xx = N(:,1:ngeom)*Xe(:,1); 
        yy = N(:,1:ngeom)*Xe(:,2); 
        zz = N*sol_e; 
        xEdge = NEdge(:,1:ngeom)*Xe(:,1); 
        yEdge = NEdge(:,1:ngeom)*Xe(:,2); 
        zEdge = NEdge*sol_e; 
        trisurf(tri,xx, yy, zz,'FaceColor','interp','EdgeColor','none')
        plot3(xEdge,yEdge,zEdge,'k'); 
    end
end
grid on
view(3); axis tight
set(gca,'FontSize',12)