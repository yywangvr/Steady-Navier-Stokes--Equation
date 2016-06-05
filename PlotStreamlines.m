function PlotStreamlines(X,velo,dom)

% Define a grid
xGrid = linspace(dom(1),dom(2),25);
yGrid = linspace(dom(3),dom(4),25);
% Interpolate the solution
tri = delaunay(X(:,1),X(:,2)); 
uGrid = tri2grid(X',tri',velo(:,1),xGrid,yGrid);   % x,y,f are column vectors. 
vGrid = tri2grid(X',tri',velo(:,2),xGrid,yGrid);   % x,y,f are column vectors. 
[xGrid,yGrid] = meshgrid(xGrid,yGrid); 

% Points where the streamlines start
y1 = min(X(:,2)); 
aux = find(abs(X(:,2)-y1) < 1e-6); 
xx = X(aux(round(length(aux)/2)),1); 
ind1 = find(abs(X(:,1)-xx) < 1e-6); 
xx = X(aux(2),1); 
ind2 = find(abs(X(:,1)-xx) < 1e-6); 
xx = X(aux(end-1),1); 
ind3 = find(abs(X(:,1)-xx) < 1e-6); 
sx = [X(ind1,1); X(ind2,1); X(ind3,1)]; 
sy = [X(ind1,2); X(ind2,2); X(ind3,2)];


h=figure(4);
if ishandle(h)
   clf(4,'reset'); 
end
figure(4); 
streamline(xGrid,yGrid,uGrid,vGrid,sx,sy)
hold on 
%plot(sx,sy,'r*')
plot(dom([1,2,2,1,1]),dom([3,3,4,4,3]),'k')
axis equal; axis tight