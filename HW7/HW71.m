close all;
clc;
file = input('please input image name: ', 's');
load(file);

source_img_coord_x = load(file).xc1;
source_img_coord_y = load(file).yc1;
target_img_coord_x = load(file).xc2;
target_img_coord_y = load(file).yc2;


xc1 = [25,50,75,80,92,75,25,30,110];
yc1 = [73,90,130,120,60,70,80,90,110];
source_img_point_num = length(source_img_coord_x);
target_img_point_num = length(target_img_coord_x);

colormap(gray);
figure(1);
subplot(2,2,1);
imagesc(source);
hold on;
plot(xc1,yc1,'*r');
hold off;
subplot(2,2,2);
imagesc(target);
hold on;
plot(xc2,yc2,'*g');
hold off;

if(method == 1)
    vector_1 = ones(source_img_point_num,1);
    vector_2 = xc1';
    vector_3 = yc1';
    vector_4 = (vector_2.*vector_3);
    vector_5 = vector_2.^2;
    vector_6 = vector_3.^2;

    matrix_A = cat(2,vector_1,vector_2,vector_3,vector_4,vector_5,vector_6);

    Wx = (matrix_A'*matrix_A)\matrix_A'*target_img_coord_x'; % polynominal parameter 
    Wy = (matrix_A'*matrix_A)\matrix_A'*target_img_coord_y';

    inserted_target_img = zeros(256,256);
        
    grid_pic = zeros(256,256);
    grid_pic(1:5:256, :)=255;
    grid_pic(:, 1:5:256)=255;
    
    inserted_grid_pic = zeros(256,256);
    


    for i = 1:256
        for j = 1:256

            target_img_coord_processed_x = Wx(1)+Wx(2)*j+Wx(3)*i+Wx(4)*j*i+Wx(5)*j^2+Wx(6)*i^2;
            target_img_coord_processed_y = Wy(1)+Wy(2)*j+Wy(3)*i+Wy(4)*j*i+Wy(5)*j^2+Wy(6)*i^2;
            if(target_img_coord_processed_x < 1 || target_img_coord_processed_y < 1 || target_img_coord_processed_x > 256 || target_img_coord_processed_y > 256)
                continue
            end

            sourse_img_gray_value = source(i,j);
            grid_pic_gray_value = grid_pic(i,j);
            inserted_target_img = insert_value(target_img_coord_processed_x, target_img_coord_processed_y,sourse_img_gray_value,inserted_target_img);
            
            inserted_grid_pic = insert_value(target_img_coord_processed_x, target_img_coord_processed_y,grid_pic_gray_value,inserted_grid_pic);


        end
    end

    for i = 1:256
        for j = 1:256
            if(inserted_target_img(i,j) > 255)
                inserted_target_img(i,j) = 255;
            end
            if(inserted_grid_pic(i,j) > 255)
                inserted_grid_pic(i,j) = 255;
            end

        end
    end



    subplot(2,2,3);
    imagesc(inserted_target_img);
    hold on;
    plot(xc2,yc2,'*g');
    hold off;
    
    
    
    subplot(2,2,4);
    imagesc(inserted_grid_pic);

    
    
    
end

if(method == 2)
    r = zeros(source_img_point_num,source_img_point_num);
    for i = 1:source_img_point_num
        for j = 1:source_img_point_num
            r(i,j) = ((xc1(i)-xc1(j))^2 + (yc1(i)-yc1(j))^2).*log(((xc1(i)-xc1(j))^2 + (yc1(i)-yc1(j))^2)+1e-8);
        end
    end
    
    vectorA1_A12 = ones(source_img_point_num,1);
    vectorB1_B12 = xc1';
    vectorC1_C12 = yc1';
    matrixA1_C12 = cat(2, vectorA1_A12,vectorB1_B12,vectorC1_C12);
    matrixA13_C15 = zeros(3,3);
    matrixA1_C15 = cat(1,matrixA1_C12,matrixA13_C15);
    matrixD1_O15 = cat(1,r,matrixA1_C12');
    matrix_A = cat(2,matrixA1_C15,matrixD1_O15);
    target_img_coord_x_zero = cat(2,target_img_coord_x,0,0,0);
    target_img_coord_y_zero = cat(2,target_img_coord_y,0,0,0);
    Wx = (matrix_A'*matrix_A)\matrix_A'*target_img_coord_x_zero'; % polynominal parameter 
    Wy = (matrix_A'*matrix_A)\matrix_A'*target_img_coord_y_zero';
    
    
    inserted_target_img = zeros(256,256);
    r_inloop = zeros(source_img_point_num,1);  
    grid_pic = zeros(256,256);
    grid_pic(1:5:256, :)=255;
    grid_pic(:, 1:5:256)=255;
    
    inserted_grid_pic = zeros(256,256);

    for i = 1:256
        for j = 1:256
            for p = 1:source_img_point_num
                r_inloop(p) = ((j-xc1(p))^2 + (i-yc1(p))^2).*log(((j-xc1(p))^2 + (i-yc1(p))^2)+1e-8);
            end
            
            target_img_coord_processed_x = Wx(1)+Wx(2)*j+Wx(3)*i+Wx(4:4+source_img_point_num-1)'*r_inloop;
            target_img_coord_processed_y = Wy(1)+Wy(2)*j+Wy(3)*i+Wy(4:4+source_img_point_num-1)'*r_inloop;
            if(target_img_coord_processed_x < 1 || target_img_coord_processed_y < 1 || target_img_coord_processed_x > 256 || target_img_coord_processed_y > 256)
                continue
            end

            sourse_img_gray_value = source(i,j);
            grid_pic_gray_value = grid_pic(i,j);
            inserted_target_img = insert_value(target_img_coord_processed_x, target_img_coord_processed_y,sourse_img_gray_value,inserted_target_img);
            
            inserted_grid_pic = insert_value(target_img_coord_processed_x, target_img_coord_processed_y,grid_pic_gray_value,inserted_grid_pic);


        end
    end

    for i = 1:256
        for j = 1:256
            if(inserted_target_img(i,j) > 255)
                inserted_target_img(i,j) = 255;
            end
            if(inserted_grid_pic(i,j) > 255)
                inserted_grid_pic(i,j) = 255;
            end

        end
    end



    subplot(2,2,3);
    imagesc(inserted_target_img);
    hold on;
    plot(xc2,yc2,'*g');
    hold off;
    
    
    
    subplot(2,2,4);
    imagesc(inserted_grid_pic);
  
end

function output_target_img = insert_value(x, y, gray_value,inserted_target_img)
output_target_img = inserted_target_img;
 x_upperbound = ceil(x);
 x_lowerbound = floor(x);
 y_upperbound = ceil(y);
 y_lowerbound = floor(y);

 diff_x_upperbound = x_upperbound - x;
 diff_y_upperbound = y_upperbound - y;

 output_target_img(y_upperbound,x_upperbound) = inserted_target_img(y_upperbound,x_upperbound) + gray_value*(1-diff_x_upperbound) * (1-diff_y_upperbound);
 output_target_img(y_upperbound,x_lowerbound) = inserted_target_img(y_upperbound,x_lowerbound) + gray_value*diff_x_upperbound * (1-diff_y_upperbound);
 output_target_img(y_lowerbound,x_upperbound) = inserted_target_img(y_lowerbound,x_upperbound) + gray_value*(1-diff_x_upperbound) * diff_y_upperbound;
 output_target_img(y_lowerbound,x_lowerbound) = inserted_target_img(y_lowerbound,x_lowerbound) + gray_value*diff_x_upperbound * diff_y_upperbound;    
end


