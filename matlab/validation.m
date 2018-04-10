%% FEM Validation %%

function validation(a, b, model, mesh, V_tm, V_te, val_flag, err_flag)

modefunc_tm = @(m, n, x, y) sin(m * pi * x / a) .* sin(n * pi * y / b);
modefunc_te = @(m, n, x, y) cos(m * pi * x / a) .* cos(n * pi * y / b);

% tm_mode_arr=[1 1; 2 1; 3 1; 1 2; 2 2];
% te_mode_arr=[1 0; 0 1; 2 0; 1 1; 2 1];


%title('Analytic Solution for Rectangular Waveguide');

%Calculate the first few modes for TE and TM cases with specified a and b
for i=1:6
    for j=1:6
        fc_mode_te_2d_arr(i,j)=((j-1)/a)^2+((i-1)/b)^2;
    end
end

fc_mode_te_1d_arr=reshape(fc_mode_te_2d_arr,[1,36]);
[fc_mode_te_1d_arr_sort, te_mode_index]=sort(fc_mode_te_1d_arr);

fc_mode_tm_2d_arr=fc_mode_te_2d_arr(2:6,2:6);
fc_mode_tm_1d_arr=reshape(fc_mode_tm_2d_arr,[1,25]);
[fc_mode_tm_1d_arr_sort, tm_mode_index]=sort(fc_mode_tm_1d_arr);

tm_mode_arr(:,1)=ceil(tm_mode_index/5);
tm_mode_arr(:,2)=mod(tm_mode_index-1,5)+1;
te_mode_arr(:,1)=floor((te_mode_index-1)/6);
te_mode_arr(:,2)=mod(te_mode_index-1,6);
te_mode_arr=te_mode_arr(2:end,:);

%In the case of a/b=2, the first five TM modes are
%TM mode1 1,1
%TM mode2 2,1
%TM mode3 3,1
%TM mode4 1,2
%TM mode5 2,2

if (val_flag==1)
    figure;
    set(gcf,'Position',[200 200 1000 500]);
end

for i = 1 : 5
    
    modeltrue =  modefunc_tm(tm_mode_arr(i,1), tm_mode_arr(i,2), mesh.Nodes(1,:), mesh.Nodes(2, :));
    %k = ((tm_mode_arr(i,1)*pi/a)^2 + (tm_mode_arr(i,2)*pi/b)^2) ^ 0.5;
    modeltrue_tm(:,i)=modeltrue / max(abs(modeltrue(:)));
    
    if (val_flag==1)
        subplot(2,5,i);
        pdeplot(model,'XYData', modeltrue_tm(:,i)')
        
        colormap jet
        xlabel('x')
        ylabel('y')
        title(['TM_{', num2str(tm_mode_arr(i,1)),num2str(tm_mode_arr(i,2)),'}']);
    end
    
end

%In the case of a/b=2, the first five TE modes are
%TE mode1 1,0
%TE mode2 0,1
%TE mode3 2,0
%TE mode4 1,1
%TE mode5 2,1
for i = 1 : 5
    
    modeltrue=modefunc_te(te_mode_arr(i,1), te_mode_arr(i,2), mesh.Nodes(1,:), mesh.Nodes(2, :));
    %k = ((te_mode_arr(i,1)*pi/a)^2 + (te_mode_arr(i,2)*pi/b)^2) ^ 0.5;
    modeltrue_te(:,i)=(modeltrue - min(modeltrue)) / (max(modeltrue) - min(modeltrue)) * 2 - 1;
    
    if (val_flag==1)
        subplot(2,5,i+5);
        pdeplot(model,'XYData', modeltrue_te(:,i)')
        
        colormap jet
        xlabel('x')
        ylabel('y')
        title(['TE_{', num2str(te_mode_arr(i,1)),num2str(te_mode_arr(i,2)),'}']);
    end
end

if (err_flag==1)
    figure;
    set(gcf,'Position',[200 200 1000 500]);
    
    for i = 1 : 5
        
        subplot(2,5,i);
        
        modelerr=modeltrue_tm - V_tm(:,i);
        meanerr=mean(modelerr);
        stderr=std(modelerr);
        
        pdeplot(model,'XYData', modelerr')
        
        colormap jet
        xlabel('x')
        ylabel('y')
        title(['TM_{', num2str(tm_mode_arr(i,1)),num2str(tm_mode_arr(i,2)),'}']);
        
    end
    
    for i = 1 : 5
        
        subplot(2,5,i+5);
        
        modelerr=modeltrue_te - V_te(:,i);
        meanerr=mean(modelerr);
        stderr=std(modelerr);
        
        pdeplot(model,'XYData', modelerr')
        
        colormap jet
        xlabel('x')
        ylabel('y')
        title(['TE_{', num2str(te_mode_arr(i,1)),num2str(te_mode_arr(i,2)),'}']);
        
    end
    
end
