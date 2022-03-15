% Main AS9
clc;
close all;

% Part1 hand registration and mean shape
figure('Name','Figure1','NumberTitle','off');
name = 'h1.mat';
subplot(2,2,1);
hold on;

num = 64;
matrix_ten_hand_x = zeros(1,num,10);
matrix_ten_hand_y = zeros(1,num,10);

for i = 1:10
    load(name);
    
    if i==10
        name='h10.mat';
    else
        name(2) = int2str(i);
    end

    image_xc1=load(name).x;
    image_yc1=load(name).y;
    plot(image_xc1, image_yc1, 'blue');
    plot(image_xc1, image_yc1, 'blue+');
    matrix_ten_hand_x(1,:,i)=image_xc1;
    matrix_ten_hand_y(1,:,i)=image_yc1;
end
hold off
title('Original shapes');

matrix_nine_realigned_hand_x=zeros(1,num,9);
matrix_nine_realigned_hand_y=zeros(1,num,9);

% 测试 1：23
% mean_nine_hand_x_mean = zeros(1,num);
% mean_nine_hand_y_mean = zeros(1,num);

% 初始化等于第一只手 测试 1：23
mean_nine_hand_x_mean = matrix_ten_hand_x(1,:,1);
mean_nine_hand_y_mean = matrix_ten_hand_y(1,:,1);

mean_ten_hand_x = zeros(1, num);
mean_ten_hand_y = zeros(1, num);

temp_mean_ten_hand_x = zeros(1,num);
temp_mean_ten_hand_y = zeros(1,num);

for j = 1:10
    
    for i=2:10
        [nine_realigned_hand_x, nine_realigned_hand_y, ~, ~] = ICP_me(matrix_ten_hand_x(1,:,i),matrix_ten_hand_y(1,:,i),mean_nine_hand_x_mean,mean_nine_hand_y_mean);% 修改matrix_ten_hand_y(1,:,1)
        matrix_nine_realigned_hand_x(1,:,i-1)=nine_realigned_hand_x;
        matrix_nine_realigned_hand_y(1,:,i-1)=nine_realigned_hand_y;
    end
    
    mean_nine_hand_x_mean = mean(matrix_nine_realigned_hand_x, 3);
    mean_nine_hand_y_mean = mean(matrix_nine_realigned_hand_y, 3);
    
    sum_nine_hand_x = sum(matrix_nine_realigned_hand_x, 3);
    sum_nine_hand_y = sum(matrix_nine_realigned_hand_y, 3);
    
    mean_ten_hand_x = (sum_nine_hand_x + matrix_ten_hand_x(1,:,1))/10;
    mean_ten_hand_y = (sum_nine_hand_y + matrix_ten_hand_y(1,:,1))/10;
    
    if (mean_ten_hand_x - temp_mean_ten_hand_x <0.1)
        break
    end
    temp_mean_ten_hand_x = mean_ten_hand_x;
    temp_mean_ten_hand_y = mean_ten_hand_y;
    
end
subplot(2,2,2);
hold on;
for i = 1:10
    if i<10
        x = matrix_nine_realigned_hand_x(1,:,i);
        y = matrix_nine_realigned_hand_y(1,:,i);
        plot(x,y, '-*b');
    end
    if i ==10
        x = matrix_ten_hand_x(1,:,1);
        y = matrix_ten_hand_y(1,:,1);
        plot(x,y, '-*b');
    end
end
hold off
% printing part
% imagesc(image_xy2);colormap(gray);
title('Realigned shapes');
subplot(2,2,3)
plot(mean_ten_hand_x, mean_ten_hand_y, 'r');
title('Mean shape');
% mean shape 不一样

% Part2 Eigenshape
mean_shape_x = zeros(1, num, 10);
mean_shape_y = zeros(1, num, 10);
for i =1:10
    mean_shape_x(1, :, i) = matrix_ten_hand_x(1,:,i) - mean_ten_hand_x;
    mean_shape_y(1, :, i) = matrix_ten_hand_y(1,:,i) - mean_ten_hand_y;
end

% put x coordinate and y coordinate together
temp_matrix_mixed_xy = cat(1, mean_ten_hand_x, mean_ten_hand_y);
temp_matrix_mixed_xy = temp_matrix_mixed_xy';

combined_de_mean_ten_hand = cat(1, mean_shape_x, mean_shape_y);

combined_de_mean_ten_hand = permute(combined_de_mean_ten_hand, [2,1,3]);

matrix_mixed_xy = reshape(combined_de_mean_ten_hand, 2*num, 10);
cov_matrix = cov(matrix_mixed_xy');

[eigenvector_matrix, ~, vector_explained_value] = pcacov(cov_matrix);
vector_explained_value = vector_explained_value ./ sum(vector_explained_value);
index_lower95 = find(cumsum(vector_explained_value) <=0.95);
index = max(index_lower95) + 1;
eigenshape = eigenvector_matrix(:,1:index);

figure('Name','Figure2','NumberTitle','off');
for i = 1:4
    subplot(2, 2, i);
    for factor = -100:10:100
        shape_disp = temp_matrix_mixed_xy + reshape(factor * eigenvector_matrix(:, i), num, 2);
        plot(shape_disp(:, 1), shape_disp(:, 2), 'b-');
        xlim([0 450]);
        ylim([0 450]);
        title(sprintf('eigen shape %d', i));
        pause(0.1);
    end
end
% subplot(2, 2, 1);
% for factor = -100:10:100
%     shape_disp = temp_matrix_mixed_xy + reshape(factor * eigenvector_matrix(:, 1), num, 2);
%     plot(shape_disp(:, 1), shape_disp(:, 2), 'b-');
%     xlim([0 450]);
%     ylim([0 450]);
%     title('First Eigenshape');
%     pause(0.1);
% end
% subplot(2, 2, 2);
% for factor = -100:10:100
%     shape_disp = temp_matrix_mixed_xy + reshape(factor * eigenvector_matrix(:, 2), num, 2);
%     plot(shape_disp(:, 1), shape_disp(:, 2), 'b-');
%     xlim([0 450]);
%     ylim([0 450]);
%     title('First Eigenshape');
%     pause(0.1);
% end
% subplot(2, 2, 3);
% for factor = -100:10:100
%     shape_disp = temp_matrix_mixed_xy + reshape(factor * eigenvector_matrix(:, 3), num, 2);
%     plot(shape_disp(:, 1), shape_disp(:, 2), 'b-');
%     xlim([0 450]);
%     ylim([0 450]);
%     title('First Eigenshape');
%     pause(0.1);
% end
% subplot(2, 2, 4);
% for factor = -100:10:100
%     shape_disp = temp_matrix_mixed_xy + reshape(factor * eigenvector_matrix(:, 4), num, 2);
%     plot(shape_disp(:, 1), shape_disp(:, 2), 'b-');
%     xlim([0 450]);
%     ylim([0 450]);
%     title('First Eigenshape');
%     pause(0.1);
% end




% part 3
% -------------------------------------------------------------------
% file = input('Input a test file: ', 's');
% %file = 'ActiveShapeTest_H2_67_HiN.mat';
% file = 'ActiveShapeTest_H7_65_noise.mat';
file = 'TestActiveShape_2019_1.mat';
% file = 'TestActiveShape_2019_2.mat';
xy = load(file).xy;
% xy = load(file).xy2;

index_matrix = zeros(index, 1);

[target_points_num,dimension_target] = size(xy);
[source_points_num,dimension_source] = size(temp_matrix_mixed_xy);

figure('Name','Figure3','NumberTitle','off');
subplot(2,2,1);
hold on;
plot(xy(:,1),xy(:,2), '+g');
plot(temp_matrix_mixed_xy(:,1),temp_matrix_mixed_xy(:,2), 'r+');
hold off;
title('Red:mean shape Green:target points');

subplot(2, 2, 2);
% match the points
for i=1:50
    new_shape = temp_matrix_mixed_xy + reshape(eigenshape*index_matrix,source_points_num,2);
    new_shape_x = new_shape(:,1)';
    new_shape_y = new_shape(:,2)';
    target_points_x = xy(:,1)';
    target_points_y = xy(:,2)';
    
    [x_icp,y_icp,T,R] = ICP_me(new_shape_x,new_shape_y,target_points_x,target_points_y);
    plot(target_points_x, target_points_y, 'g*-');
    hold on;
    plot(x_icp, y_icp, 'b*-');
    hold off;
    xlim([0 450]);
    ylim([0 450]);
    
    [matched_point_x,matched_point_y] = match_point(new_shape_x,new_shape_y,target_points_x,target_points_y);
%     [matched_point_x, matched_point_y] = bestmatch(new_shape, xy);
    nearest_point_xy = cat(1, matched_point_x,matched_point_y)';
    
    nearest_point_xy_inverted = (R \ (nearest_point_xy - repmat(T, 1, source_points_num)')')';
    nearest_point_x_inverted = nearest_point_xy_inverted(:, 1)';
    nearest_point_y_inverted = nearest_point_xy_inverted(:, 2)';
    
    index_matrix = eigenshape' * reshape((cat(1, nearest_point_x_inverted, nearest_point_y_inverted)' - temp_matrix_mixed_xy), 2 * source_points_num, 1);
    pause(0.5);
end

xlim([0 450]);
ylim([0 450]);

new_shape = temp_matrix_mixed_xy + reshape(eigenshape * index_matrix, source_points_num, 2);
x_new_shape = new_shape(:, 1)';
y_new_shape = new_shape(:, 2)';
[x_new_shape, y_new_shape, T, R] = ICP(x_new_shape, y_new_shape, target_points_x, target_points_y);

subplot(2, 2, 3);

plot(target_points_x, target_points_y, 'g*', target_points_x, target_points_y, 'g-');
hold on;
plot(x_new_shape, y_new_shape, 'b*', x_new_shape, y_new_shape, 'b-');
hold off;
xlim([0 450]);
ylim([0 450]);
title('Final results');

