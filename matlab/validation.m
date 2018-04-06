%% FEM Validation 

function validation(a, b, model, mesh)

% if (a-b*2>0.0000001)
%     disp('a=2*b condition not satisfied!');
%     return;
% end

modefunc_tm = @(m, n, x, y) sin(m * pi * x / a) .* sin(n * pi * y / b);
modefunc_te = @(m, n, x, y) cos(m * pi * x / a) .* cos(n * pi * y / b);

tm_mode_arr=[1 1; 2 1; 3 1; 1 2; 2 2];
te_mode_arr=[1 0; 0 1; 2 0; 1 1; 2 1];

figure;
set(gcf,'Position',[200 200 1200 600]);
%title('Analytic Solution for Rectangular Waveguide');
%TM mode1 1,1
%TM mode2 2,1 
%TM mode3 3,1
%TM mode4 1,2
%TM mode5 2,2
for i = 1 : 5
    
    subplot(2,5,i);
    modeltrue =  modefunc_tm(tm_mode_arr(i,1), tm_mode_arr(i,2), mesh.Nodes(1,:), mesh.Nodes(2, :));
    modeltrue=(modeltrue - min(modeltrue)) / (max(modeltrue) - min(modeltrue)) * 2 - 1;
    pdeplot(model,'XYData', modeltrue')

    colormap jet
    xlabel('x')
    ylabel('y')
    title(['TM_{', num2str(tm_mode_arr(i,1)),num2str(tm_mode_arr(i,2)),'}']);

end

%TE mode1 1,0
%TE mode2 0,1 
%TE mode3 2,0
%TE mode4 1,1
%TE mode5 2,1
for i = 1 : 5
    
    subplot(2,5,i+5);
    modeltrue =  modefunc_te(te_mode_arr(i,1), te_mode_arr(i,2), mesh.Nodes(1,:), mesh.Nodes(2, :));
    modeltrue=(modeltrue - min(modeltrue)) / (max(modeltrue) - min(modeltrue)) * 2 - 1;
    pdeplot(model,'XYData', modeltrue')

    colormap jet
    xlabel('x')
    ylabel('y')
    title(['TE_{', num2str(te_mode_arr(i,1)),num2str(te_mode_arr(i,2)),'}']);

end