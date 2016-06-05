function [X,T] = CreateUniformMesh(dom,nx,ny,elem,degree)
% [X,T] = CreateUniformMesh(dom,nx,ny,elem,degree)
% Uniform mesh in a rectangular domain
% Input:
%   dom = [x1,x2,y1,y2]:  vertices' coordinates
%   nx,ny: number of elements in each direction
%   elem: type of element (0:quadrilateral, 1:triangle, 11: triangle with bubble function)
%   degree: interpolation degree
% Output:
%   X:  nodal coordinates
%   T:  connectivities

x1 = dom(1); x2 = dom(2);
y1 = dom(3); y2 = dom(4);

npx = degree*nx + 1;
npy = degree*ny + 1;

npt = npx*npy;
x = linspace(x1,x2,npx);
y = linspace(y1,y2,npy);
[x,y] = meshgrid(x,y);
X = [reshape(x',npt,1), reshape(y',npt,1)];

if elem == 0
    nen = (degree+1)^2;
    T = zeros(nx*ny,nen);
    if degree == 1
        for i=1:ny
            for j=1:nx
                ielem = (i-1)*nx+j;
                inode = (i-1)*(npx)+j;
                T(ielem,:) = [inode   inode+1   inode+npx+1   inode+npx];
            end
        end
    elseif degree == 2
        for i=1:ny
            for j=1:nx
                ielem = (i-1)*nx + j;
                inode = (i-1)*2*npx + 2*(j-1) + 1;
                nodes_aux = [inode+(0:2)  inode+npx+(0:2)  inode+2*npx+(0:2)];
                T(ielem,:) = nodes_aux([1  3  9  7  2  6  8  4  5]);
            end
        end
    else
        error('not available element')
    end
elseif elem == 1
    nen = (degree+1)*(degree+2)/2;
    T = zeros(2*nx*ny,nen);
    if degree == 1
        for i=1:ny
            for j=1:nx
                ielem = 2*((i-1)*nx+j)-1;
                inode = (i-1)*(npx)+j;
                T(ielem,:) = [inode   inode+1   inode+(npx)];
                T(ielem+1,:) = [inode+1   inode+1+npx   inode+npx];
            end
        end
        % Modification of left lower and right upper corner elements to avoid them
        % having all their nodes on the boundary
        if npx > 2
            T(1,:) = [1  npx+2   npx+1];
            T(2,:) = [1    2     npx+2];
            aux = size(T,1);
            T(aux-1,:) = [npx*ny-1    npx*npy   npx*npy-1];
            T(aux,:)   = [npx*ny-1    npx*ny    npx*npy];
        end
    elseif degree == 2
        for i=1:ny
            for j=1:nx
                ielem=2*((i-1)*nx+j)-1;
                inode=(i-1)*2*(npx)+2*(j-1)+1;
                nodes_aux = [inode+(0:2)  inode+npx+(0:2)  inode+2*npx+(0:2)];
                T(ielem,:) = nodes_aux([1  3  7  2  5  4]);
                T(ielem+1,:) = nodes_aux([3  9  7  6  8  5]);
            end
        end
        % Modification of left lower and right upper corner elements to avoid them
        % having all their nodes on the boundary
        if npx > 3
            inode = 1;
            nodes_aux = [inode+(0:2)  inode+npx+(0:2)  inode+2*npx+(0:2)];
            T(1,:) = nodes_aux([1  9  7  5  8  4]);
            T(2,:) = nodes_aux([1  3  9  2  6  5]);
            
            ielem  = size(T,1)-1;
            inode = npx*(npy-2)-2;
            nodes_aux = [inode+(0:2)  inode+npx+(0:2)  inode+2*npx+(0:2)];
            T(ielem,:) = nodes_aux([1  9  7  5  8  4]);
            T(ielem+1,:) = nodes_aux([1  3  9  2  6  5]);
        end
    else
        error('not available element')
    end
elseif elem == 11
    if degree == 1
        T = zeros(2*nx*ny,4);
        npt = size(X,1); 
        for i=1:ny
            for j=1:nx
                ielem = 2*((i-1)*nx+j)-1;
                inode = (i-1)*(npx)+j;
                n_ad = npt + 2*((i-1)*nx+j)-1;
                T(ielem,:) = [inode   inode+1   inode+(npx)   n_ad];
                T(ielem+1,:) = [inode+1   inode+1+npx   inode+npx   n_ad+1];
            end
        end
        % Modification of left lower and right upper corner elements to avoid them
        % having all their nodes on the boundary
        if npx > 2
            T(1,:) = [1  npx+2   npx+1   npt+1];
            T(2,:) = [1    2     npx+2   npt+2];
            aux = size(T,1);
            T(aux-1,:) = [npx*ny-1    npx*npy   npx*npy-1   npt+aux-1];
            T(aux,:)   = [npx*ny-1    npx*ny    npx*npy     npt+aux];
        end
    else
        error('not available element')
    end
else
    error('not available element')
end




% elseif elem == 1
%     nen = (degree+1)*(degree+2)/2;
%     T = zeros(2*nx*ny,nen);
%     nx_2 = round(nx/2); ny_2 = round(ny/2);
%     if degree == 1
%         for i=1:ny
%             for j=1:nx
%                 ielem = 2*((i-1)*nx+j)-1;
%                 inode = (i-1)*(npx)+j;
%                 nodes = [inode   inode+1   inode+npx+1    inode+npx];
%                 if (i<=ny_2 && j<=nx_2) || (i>ny_2 && j>nx_2)
%                     T(ielem,:) = nodes([1,2,3]);
%                     T(ielem+1,:) = nodes([1,3,4]);
%                 else
%                     T(ielem,:) = nodes([1,2,4]);
%                     T(ielem+1,:) = nodes([2,3,4]);
%                 end
%             end
%         end
%     elseif degree == 2
%         for i=1:ny
%             for j=1:nx
%                 ielem=2*((i-1)*nx+j)-1;
%                 inode=(i-1)*2*(npx)+2*(j-1)+1;
%                 nodes = [inode+(0:2)  inode+npx+(0:2)  inode+2*npx+(0:2)];
%                 if (i<=ny_2 && j<=nx_2) || (i>ny_2 && j>nx_2)
%                     T(ielem,:)   = nodes([1  3  9  2  6  5]);
%                     T(ielem+1,:) = nodes([1  9  7  5  8  4]);
%                 else
%                     T(ielem,:)   = nodes([1  3  7  2  5  4]);
%                     T(ielem+1,:) = nodes([3  9  7  6  8  5]);
%                 end
%             end
%         end
%
%     else
%         error('not available element')
%     end
% elseif elem == 11
%     if degree == 1
%         T = zeros(2*nx*ny,4);
%         nx_2 = round(nx/2); ny_2 = round(ny/2);
%         for i=1:ny
%             for j=1:nx
%                 ielem = 2*((i-1)*nx+j)-1;
%                 inode = (i-1)*(npx)+j;
%                 nodes = [inode   inode+1   inode+npx+1    inode+npx];
%                 n_ad = npx*npy + 2*((i-1)*nx+j)-1;
%                 if (i<=ny_2 && j<=nx_2) || (i>ny_2 && j>nx_2)
%                     T(ielem,:) = [nodes([1,2,3]), n_ad];
%                     T(ielem+1,:) = [nodes([1,3,4]), n_ad+1];
%                 else
%                     T(ielem,:) = [nodes([1,2,4]), n_ad];
%                     T(ielem+1,:) = [nodes([2,3,4]), n_ad+1];
%                 end
%             end
%         end
%     else
%         error('not available element')
%     end