clear
clc

minx = 0;
miny = 0;
a = 1;
b = 0.5;
maxx = a;
maxy = b;

ns = [82;49;0];
gd = [3;4;minx;minx;maxx;maxx;miny;maxy;maxy;miny];
sf = 'R1';
dl = decsg(gd, sf, ns);

model = createpde(1);
geometryFromEdges(model, dl);
mesh = generateMesh(model, 'Hmax', 0.02);

xdata = mesh.Nodes(1, :);
ydata = mesh.Nodes(2, :);
tol = 1e-8 * a;

dirichlet = abs(xdata - 0) < tol | abs(xdata - a) < tol | abs(ydata - 0) < tol | abs(ydata - b) < tol;
% dirichlet =  abs(xdata - 0) < tol & abs(ydata - 0) < tol;
figure(1)

pdeplot(model)
hold on
plot(mesh.Nodes(1,dirichlet),mesh.Nodes(2,dirichlet),'or','MarkerFaceColor','g')
xlabel('x')
ylabel('y')
axis equal

% prompt = 'Mesh OK?';
% x = input(prompt);
% if isempty(x)

%%
disp('Generating Mesh File');
fileID = fopen('../output/mesh.bin', 'wb');
if fileID == -1
    cd ..
    !mkdir output
    cd matlab
    fileID = fopen('../output/mesh.bin', 'wb');
end

numNode = size(mesh.Nodes, 2);
numEm = size(mesh.Elements, 2);

fwrite(fileID, size(mesh.Nodes, 2), 'int');
fwrite(fileID, mesh.Nodes, 'double');
fwrite(fileID, size(mesh.Elements, 2), 'int');
fwrite(fileID, mesh.Elements, 'int');
fclose(fileID);

cd ..
!ece618.exe
cd matlab
%%
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
num_eigval = 4;
[eigfunc, eigval] = eigs(A(~dirichlet, ~dirichlet), B(~dirichlet, ~dirichlet), num_eigval, 'sm');
V = zeros([numNode num_eigval]);
V(~dirichlet, :) = eigfunc;


V = fliplr(V);
eigval = fliplr(eigval);
modefunc = @(m, n, x, y) sin(m * pi * x / a) .* sin(n * pi * y / b);
modenum = 4;
mode = V(:, modenum);
modetrue =  modefunc(1, 2, mesh.Nodes(1,:), mesh.Nodes(2, :));

figure(2)
mode = mode / max(abs(mode(:)));
pdeplot(model,'ZData', mode)
xlabel('x')
ylabel('y')
title(['Mode ', num2str(modenum)])


%% TEz
num_eigval = 6;
[eigfunc, eigval] = eigs(A, B, num_eigval, 'sm');

V = fliplr(eigfunc);
eigval = fliplr(eigval);
modefunc = @(m, n, x, y) cos(m * pi * x / a) .* cos(n * pi * y / b);
modenum = 5;
mode = V(:, modenum + 1);
mode1true =  modefunc(1, 1, mesh.Nodes(1,:), mesh.Nodes(2, :));

figure(3)
%normalize the functions
span = max(mode(:)) - min(mode(:));
mode = (mode - min(mode(:))) / span * 2 - 1;

pdeplot(model,'ZData', mode)
xlabel('x')
ylabel('y')
title(['Mode ', num2str(modenum)])
