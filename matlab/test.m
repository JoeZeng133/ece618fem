% clear
% clc

a = 1; %width
b = 0.5; %height
order = 1;
em_size = 0.1;

minx = 0;
miny = 0;
maxx = a;
maxy = b;

% generation of a rectangle 
ns = [82;49;0];
gd = [3;4;minx;minx;maxx;maxx;miny;maxy;maxy;miny];
sf = 'R1';
dl = decsg(gd, sf, ns);

model = createpde(1);
geometryFromEdges(model, dl);
% 'Hmax' is the mesh size, 'GeomericOrder' is the mesh element order
% for 2nd order mesh, let 'GeometricOrder'='quadratic' and call quad.exe
% for 1st order mesh, let 'GeometricOrder'='linear'    and call linear.exe
if order == 1
    mesh = generateMesh(model, 'Hmax', 0.1, 'GeometricOrder', 'linear'); %generate mesh with maximum element length 0.02
else
    mesh = generateMesh(model, 'Hmax', 0.1, 'GeometricOrder', 'quadratic'); %generate mesh with maximum element length 0.02
end


xdata = mesh.Nodes(1, :);
ydata = mesh.Nodes(2, :);
tol = 1e-8 * a;

%look for boundary nodes
dirichlet = abs(xdata - 0) < tol | abs(xdata - a) < tol | abs(ydata - 0) < tol | abs(ydata - b) < tol;

figure(1)
pdeplot(model), hold on
plot(mesh.Nodes(1,dirichlet),mesh.Nodes(2,dirichlet),'or','MarkerFaceColor','g'), hold off
xlabel('x')
ylabel('y')
axis equal


%% generate mesh.bin for assembler
disp('Generating Mesh File');
fileID = fopen('../output/mesh.bin', 'wb');
if fileID == -1
    cd ..
    !mkdir output
    cd matlab
    fileID = fopen('../output/mesh.bin', 'wb');
end

numNode = size(mesh.Nodes, 2); %number of nodes
numEm = size(mesh.Elements, 2); %number of elements

fwrite(fileID, size(mesh.Nodes, 2), 'int'); 
fwrite(fileID, mesh.Nodes, 'double'); %nodes coordinate 
fwrite(fileID, size(mesh.Elements, 2), 'int'); 
fwrite(fileID, mesh.Elements, 'int'); %element connection
fclose(fileID);

% run assembler, for 2nd order mesh, call quad.exe instead
cd ..
if order == 1
    !linear.exe 
else
    !quad.exe
end

cd matlab
%% read sparse matrix A and B
fileID = fopen('../output/equation.bin', 'rb');
if fileID == -1
    disp('Programm yet to be run');
end

numEqs = fread(fileID, 1, 'int');
Ai = fread(fileID, numEqs, 'int');
Aj = fread(fileID, numEqs, 'int');
Av = fread(fileID, numEqs, 'double');

Bi = fread(fileID, numEqs, 'int');
Bj = fread(fileID, numEqs, 'int');
Bv = fread(fileID, numEqs, 'double');


A = sparse(Ai, Aj, Av, numNode, numNode);
B = sparse(Bi, Bj, Bv, numNode, numNode);
fclose(fileID);

%% TMz
% solve generalized eigenvalue problem (A - lambda * B) * x = 0
num_eig = 6;
Aminusdir = A(~dirichlet, ~dirichlet);
Bminusdir = B(~dirichlet, ~dirichlet);

[eigfunc, eigval] = eigs(Aminusdir, Bminusdir, num_eig , 'SA'); %impose dirichlet boundary
eigval = diag(eigval);

V = zeros([numNode num_eig]);
V(~dirichlet, :) = eigfunc;

% modefunc = @(m, n, x, y) sin(m * pi * x / a) .* sin(n * pi * y / b);
% mode5 =  modefunc(4, 1, mesh.Nodes(1,:), mesh.Nodes(2, :))';
% mode6 =  modefunc(2, 2, mesh.Nodes(1,:), mesh.Nodes(2, :))';
 
figure(2)
title('TM Mode')
for i = 1 : num_eig
    subplot(2, 3,i)
    mode = V(:, i);
    mode = mode / max(abs(mode(:)));
    pdeplot(model,'XYData', mode)
    colormap jet
    xlabel('x')
    ylabel('y')
    axis equal
    axis([0 a 0 b])
    title(['TM k= ', num2str(eigval(i))])
end


%% TEz
num_eig = 6;
[eigfunc, eigval, flag] = eigs(A,B,num_eig + 1,'sa');

V = eigfunc;
eigval = diag(eigval);

% modefunc = @(m, n, x, y) cos(m * pi * x / a) .* cos(n * pi * y / b);
% mode2 =  modefunc(2, 0, mesh.Nodes(1,:), mesh.Nodes(2, :))';
% mode3 =  modefunc(0, 1, mesh.Nodes(1,:), mesh.Nodes(2, :))';

figure(3)
title('TE Mode')
for i = 1 : num_eig
    subplot(2, 3,i)
    mode = V(:, i + 1);
    mode = (mode - min(mode)) / (max(mode) - min(mode)) * 2 - 1;
    pdeplot(model,'XYData', mode)
    colormap jet
    xlabel('x')
    ylabel('y')
    axis equal
    axis([0 a 0 b])
    title(['TE k =  ', num2str(eigval(i + 1))])
end