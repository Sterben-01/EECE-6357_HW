close all;
clc;
figure(1);
subplot(2,2,1);
hold on;

h_x = zeros(1,64,10);%raw data 10 hand
h_y = zeros(1,64,10);



for i = 1:10
    ii = int2str(i);
    h = 'hh';
    mat = '.mat';
    file = append(h, ii, mat);
    load(file);
    plot(x,y,'-*b');
    h_x(1,:,i) = x;
    h_y(1,:,i) = y;
end


hold off;

registured_hand_x_list = zeros(1,64,9);
registured_hand_y_list = zeros(1,64,9);

mean_9_hand_x = zeros(1,64);
mean_9_hand_y = zeros(1,64);

first_hand_x = h_x(1,:,1);
first_hand_y = h_x(1,:,1);

subplot(2,2,2);

counter = 0;


pre_mean_10_hand_x = zeros(1,64);
pre_mean_10_hand_y = zeros(1,64);

for oo = 1:100
    
    for i = 2:10
        [registured_hand_x, registured_hand_y, ~, ~] = ICPme(h_x(1,:,i), h_y(1,:,i),first_hand_x ,first_hand_y);
%         [registured_hand_x, registured_hand_y, ~, ~] = ICP(h_x(1,:,i), h_y(1,:,i),h_x(1,:,1), h_y(1,:,1));
        registured_hand_x_list(1,:,i-1) = registured_hand_x;
        registured_hand_y_list(1,:,i-1) = registured_hand_y;
        
    end
    

    mean_9_hand_x = mean(registured_hand_x_list,3);
    mean_9_hand_y = mean(registured_hand_y_list,3);
    
    
    first_hand_x = mean_9_hand_x;
    first_hand_y = mean_9_hand_y;
    
    sum_9_hand_x = sum(registured_hand_x_list,3);
    sum_9_hand_y = sum(registured_hand_y_list,3);

    
    mean_10_hand_x = (sum_9_hand_x+ h_x(1,:,1)) / 10;
    mean_10_hand_y = (sum_9_hand_y+ h_y(1,:,1)) / 10;
    
    
    if(norm(mean_10_hand_x - pre_mean_10_hand_x) < 0.1)
        
        break
    end
    
    pre_mean_10_hand_x = mean_10_hand_x;
    pre_mean_10_hand_y = mean_10_hand_y;
    counter = counter +1;
    
end

hold on;
for i = 1:9
    x = registured_hand_x_list(1,:,i);
    y = registured_hand_y_list(1,:,i);
    plot(x,y,'-*b');
end
hold off;


subplot(2,2,3);
plot(mean_10_hand_x, mean_10_hand_y, 'r');

%phase2

de_mean_10_hand_x = zeros(1,64,10); 
de_mean_10_hand_y = zeros(1,64,10); 

for i = 1:10
    de_mean_10_hand_x(1,:,i) = h_x(1,:,i) - mean_10_hand_x;

    de_mean_10_hand_y(1,:,i) = h_y(1,:,i) - mean_10_hand_y;
end

combined_mean_10_hand = cat(1,mean_10_hand_x, mean_10_hand_y);

combined_mean_10_hand = combined_mean_10_hand';

combined_de_mean_10_hand = cat(1,de_mean_10_hand_x, de_mean_10_hand_y);


combined_de_mean_10_hand = permute(combined_de_mean_10_hand,[2,1,3]);



reshape_combined_de_mean_10_hand = reshape(combined_de_mean_10_hand,128,10);

covariance_matrix = cov(reshape_combined_de_mean_10_hand');

[column, ~, exp] = pcacov(covariance_matrix);

exp = exp ./ sum(exp);

index_less_p95 = find(cumsum(exp) <= 0.95);

index_index_less_p95 = max(index_less_p95) + 1;

eigen = column(:, 1:index_index_less_p95);

figure(2)

subplot(2,2,1);
for i = -100:10:100
    contour = combined_mean_10_hand + reshape(i * column(:,1),64,2);
    plot(contour(:, 1), contour(:, 2), 'b-');
    pause(0.1);
    
end
subplot(2,2,2);
for i = -100:10:100
    contour = combined_mean_10_hand + reshape(i * column(:,2),64,2);
    plot(contour(:, 1), contour(:, 2), 'b-');
    pause(0.1);
    
end

subplot(2,2,3);
for i = -100:10:100
    contour = combined_mean_10_hand + reshape(i * column(:,3),64,2);
    plot(contour(:, 1), contour(:, 2), 'b-');
    pause(0.1);
    
end

subplot(2,2,4);
for i = -100:10:100
    contour = combined_mean_10_hand + reshape(i * column(:,4),64,2);
    plot(contour(:, 1), contour(:, 2), 'b-');
    pause(0.1);
    
end



%phase 3


inport_new_point = input('The name of your test file is: ', 's');
xy = load(inport_new_point).xy2;

b = zeros(index_index_less_p95, 1);

[new_points_num, ~] = size(xy);
[~,mean_shape_points_num] = size(mean_10_hand_x);

figure(3)
subplot(2,2,1);
plot(xy(:, 1), xy(:, 2), 'g*-');

hold on;

plot(mean_10_hand_x, mean_10_hand_y, 'r');

hold off;

subplot(2,2,2);

for i = 1:20
    
    contour_b = combined_mean_10_hand + reshape(eigen * b, mean_shape_points_num, 2);
    contour_x = contour_b(:, 1)';
    contour_y = contour_b(:, 2)';
    new_point_x = xy(:, 1)';
    new_point_y = xy(:, 2)';
    
    [shape_x, shape_y, T, r] = ICPme(contour_x, contour_y, new_point_x, new_point_y);
    
    plot(new_point_x, new_point_y, 'g*-');
    hold on;
    plot(shape_x, shape_y, 'b*-');
    hold off;
    
    [calced_matrix_x, calced_matrix_y] = distance_match(contour_x, contour_y, new_point_x, new_point_y);
%     [calced_matrix_x, calced_matrix_y] = bestmatch(contour_b, xy);
    
    
    combined_calced_matrix_xy = cat(1, calced_matrix_x, calced_matrix_y)';
    
    combined_calced_matrix_xy_inv = (r\(combined_calced_matrix_xy - repmat(T,1,mean_shape_points_num)')')';
    
    combined_calced_matrix_x_inv = combined_calced_matrix_xy_inv(:,1)';
    combined_calced_matrix_y_inv  = combined_calced_matrix_xy_inv(:,2)';
    
    b = eigen' * reshape((cat(1,combined_calced_matrix_x_inv, combined_calced_matrix_y_inv)' - combined_mean_10_hand), mean_shape_points_num*2, 1);  
    
    pause(0.5);
    
    
%     plot(new_point_x, new_point_y, 'g*-');
%     hold on;
%     plot(shape_x, shape_y, 'b*-');
%     plot(calced_matrix_x(1,:),calced_matrix_y(1,:),'-r*');
%     hold off;

end

contour_b = combined_mean_10_hand + reshape(eigen * b, mean_shape_points_num, 2);
contour_x = contour_b(:, 1)';
contour_y = contour_b(:, 2)';
[shape_x, shape_y, T, r] = ICPme(contour_x, contour_y, new_point_x, new_point_y);
subplot(2,2,3);
plot(new_point_x, new_point_y, 'g*-');
hold on;
plot(shape_x, shape_y, 'r');
hold off;






function [processed_point_x, processed_point_y, T, r] = ICPme (initial_point_x, initial_point_y, target_point_x, target_point_y)
    
    blue_point_number = length(initial_point_x);
    
    
    T = [0; 0];
    r = [1, 0; 0, 1];

    for iter = 1:500
        red_point_x = zeros(1,blue_point_number);
        red_point_y = zeros(1,blue_point_number);

        for i = 1:blue_point_number
            min_blue_value = 100000000;
            min_x = 100000;
            min_y = 100000;

            for j = 1:length(target_point_x)
                distance_squre = (initial_point_x(i) - target_point_x(j))^2 + (initial_point_y(i) - target_point_y(j))^2;
                if(distance_squre < min_blue_value)
                    min_blue_value = distance_squre;
                    min_x = target_point_x(j);
                    min_y = target_point_y(j);
                end        
            end

            red_point_x(i) = min_x;
            red_point_y(i) = min_y;

        end

        blue_point_x_mean = mean(initial_point_x);
        blue_point_y_mean = mean(initial_point_y);
        red_point_x_mean = mean(red_point_x);
        red_point_y_mean = mean(red_point_y);


        anti_mean_blue_point_x = initial_point_x - blue_point_x_mean;
        anti_mean_blue_point_y = initial_point_y - blue_point_y_mean;
        anti_mean_red_point_x = red_point_x - red_point_x_mean;
        anti_mean_red_point_y = red_point_y - red_point_y_mean;


        combined_blue_point = cat(1, anti_mean_blue_point_x, anti_mean_blue_point_y);
        combined_red_point = cat(1, anti_mean_red_point_x, anti_mean_red_point_y);


        H = combined_blue_point * combined_red_point';
        [U, ~, V] = svd(H);

        R = V * U';

        t = [red_point_x_mean; red_point_y_mean] - R * [blue_point_x_mean; blue_point_y_mean];

%         transformed_blue_point = R * combined_blue_point + t + [blue_point_x_mean; blue_point_y_mean];
        transformed_blue_point = t + R * cat(1, initial_point_x, initial_point_y);
        

        
        transformed_blue_point_x = transformed_blue_point(1,:);
        transformed_blue_point_y = transformed_blue_point(2,:);
        
        T = T + t;
        r = r * R;
        
        initial_point_x = transformed_blue_point_x;
        initial_point_y = transformed_blue_point_y;

    end
    processed_point_x = initial_point_x;
    processed_point_y = initial_point_y;
    
end


function [calced_matrix_x, calced_matrix_y] = distance_match(shape_x, shape_y,target_x, target_y)
    %calc zhou chang
    perimeter_shape = 0;
    perimeter_target = 0;
    
    calced_matrix_x = zeros(1, length(shape_x));
    
    calced_matrix_y = zeros(1, length(shape_y));
    
    
    for i = 1:length(shape_x)-1
        perimeter_shape = perimeter_shape + norm([(shape_x(i+1) - shape_x(i)), (shape_y(i+1) - shape_y(i))] ,2);
    end
    
    for i = 1:length(target_x)-1
        perimeter_target = perimeter_target + norm([(target_x(i+1) - target_x(i)), (target_y(i+1) - target_y(i))] ,2);
    end
    
    cumulate_shape_used_distance = 0;
    
    for j = 2:length(shape_x)-1
        
        temp_distance = norm([(shape_x(j) - shape_x(j-1)), (shape_y(j) - shape_y(j-1))] ,2);
        cumulate_shape_used_distance = cumulate_shape_used_distance + temp_distance;
       
        cumulate_target_used_distance = 0;
        
        
        
        for p = 1:length(target_x)-1
            
            temp_distance_target = norm([(target_x(p+1) - target_x(p)), (target_y(p+1) - target_y(p))] ,2);
            
            
            if((cumulate_target_used_distance + temp_distance_target) / perimeter_target) >= (cumulate_shape_used_distance / perimeter_shape)
                
                A = cumulate_shape_used_distance / perimeter_shape; 
                B = cumulate_target_used_distance / perimeter_target;
                C = (cumulate_target_used_distance + temp_distance_target) / perimeter_target;
                
                
                delta = (cumulate_target_used_distance + temp_distance_target)/perimeter_target - cumulate_shape_used_distance/perimeter_shape;
                portion = 1 - delta / (temp_distance_target / perimeter_target);
                
%                 disp(A)
%                 disp(B)
%                 disp(C)
                precentage_in_use = (A-B)/(C-B);
                disp(p)
                disp(precentage_in_use);
                
                calced_matrix_x(j) = target_x(p) + portion*(target_x(p+1) - target_x(p));
                calced_matrix_y(j) = target_y(p) + portion*(target_y(p+1) - target_y(p));
                
                break;
                
            end
            
            cumulate_target_used_distance = cumulate_target_used_distance + temp_distance_target;
            
        end

    end
    
    calced_matrix_x(1) = target_x(1);
    calced_matrix_y(1) = target_y(1);
    calced_matrix_x(length(shape_x)) = target_x(length(target_x));
    calced_matrix_y(length(shape_y)) = target_y(length(target_y));
    
end
















