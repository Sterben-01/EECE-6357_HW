close all;
clc;
file = input('please input image name: ', 's');
load(file);

static = load(file).static;
dynamic = load(file).dynamic;
smooth = load(file).std;

figure(1);
colormap(gray)
subplot(2,2,1);
imagesc(static);
subplot(2,2,2);
imagesc(dynamic);
resize_value = 256;

if(numlevel > 1 )
    for i = 1:numlevel-1
        resize_value = resize_value/2;
    end
    
end
static_resize = imresize(static,[resize_value,resize_value]);
dynamic_resize = imresize(dynamic,[resize_value,resize_value]);

dynamic_resize_after = dynamic_resize;

[static_gradient_y,static_gradient_x] = gradient(static_resize);


transformation_field_x=zeros(size(dynamic_resize)); 
transformation_field_y=zeros(size(dynamic_resize));

% smooth_filter = fspecial('gaussian',35,smooth);
disp(resize_value);
for j = 1:numlevel

    for i=1:200

            image_diff=dynamic_resize_after-static_resize;


            Vx = -(image_diff.*static_gradient_x)./((static_gradient_x.^2+static_gradient_y.^2)+image_diff.^2);
            Vy = -(image_diff.*static_gradient_y)./((static_gradient_x.^2+static_gradient_y.^2)+image_diff.^2);


            Vx(isnan(Vx))=0; 
            Vy(isnan(Vy))=0;
            smooth_filter_x = imgaussfilt(Vx, smooth, 'FilterDomain', 'spatial');
            smooth_filter_y = imgaussfilt(Vy, smooth, 'FilterDomain', 'spatial');

    %         Vx_smooth=imfilter(Vx,smooth_filter);
    %         Vy_smooth=imfilter(Vy,smooth_filter);
    % 
    %         transformation_field_x=transformation_field_x+Vx_smooth;
    %         transformation_field_y=transformation_field_y+Vy_smooth;
            transformation_field_x=transformation_field_x+smooth_filter_x;
            transformation_field_y=transformation_field_y+smooth_filter_y;

    % 
    %         meshgrid_x = 1:256;
    %         meshgrid_y = 1:256;
    %         
    %         [meshgrid_X,meshgrid_Y] = meshgrid(meshgrid_x,meshgrid_y);
    %         meshgrid_X_add = meshgrid_X + transformation_field_x;
    %         meshgrid_Y_add = meshgrid_Y + transformation_field_y;
    %         
    %         dynamic_resize = interp2(meshgrid_X, meshgrid_Y, dynamic_resize, meshgrid_X_add, meshgrid_Y_add);
            dynamic_resize_after = movepixels(dynamic_resize, transformation_field_x, transformation_field_y);
            pause(0.01);
            subplot(2,2,3);
            imagesc(dynamic_resize_after);
            pic_four = abs(static_resize - dynamic_resize_after);
            subplot(2,2,4);
            imagesc(pic_four);
            
    end
    resize_value = resize_value * 2;
    static_resize = imresize(static,[resize_value,resize_value]);
    dynamic_resize = imresize(dynamic,[resize_value,resize_value]);
    
    dynamic_resize_after = dynamic_resize;

    [static_gradient_y,static_gradient_x] = gradient(static_resize);
    transformation_field_x=imresize(transformation_field_x,[resize_value,resize_value]); 
    transformation_field_y=imresize(transformation_field_y,[resize_value,resize_value]);
    disp(resize_value);
end





