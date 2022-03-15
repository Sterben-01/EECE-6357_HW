close all;
clc;
file = input('please input image name: ', 's');
load(file);
background = zeros(256,256);
imagesc(background);
colormap summer;
figure(1);
hold on;
plot(xc, yc, 'r');
plot(xc, yc, '*r');
plot(pointx, pointy, 'g');
plot(pointx, pointy, '*b');
hold off;



% pic2 = double(roipoly(background, xc, yc));
% edge_image = edge(pic2, 'Canny');
% distance_map = double(bwdist(edge_image));


pre_edge_image = zeros(256);
[resample_x, resample_y] = resample(500, xc, yc);
figure(5);
plot(resample_x, resample_y, 'g+');
set(gca,'YDir','reverse');
edge_image = create_edge(resample_y, resample_x, pre_edge_image);
distance_map = double(bwdist(edge_image));
figure(2);
hold on;
set(gca,'YDir','reverse');
imagesc(distance_map);
plot(pointx, pointy, 'g');
plot(pointx, pointy, '*g');
hold off;


if(method == 1)
    param = [0,0,0];

    figure(3);
    hold on;
    set(gca,'YDir','reverse');
    imagesc(distance_map);
    options = optimset('MaxIter', 10000, 'TolFun', 1e-11, 'TolX', 1e-11);
    optim = fminsearch(@angle_translation, param, options);
    hold off;


    v = cat(1,pointx,pointy);
    x_center = 128;
    y_center = 128;
    rotation_center = repmat([x_center; y_center], 1, length(pointx));
    rotation_matrix = [cos(optim(1)) -sin(optim(1)); sin(optim(1)) cos(optim(1))];
    rotation_center_translation = repmat([optim(2); optim(3)], 1, length(pointx));
    shift_point = v - rotation_center;     
    shift_point_applied_origion = rotation_matrix*shift_point;           
    shift_point_back = shift_point_applied_origion + rotation_center + rotation_center_translation;  

    rotated_x = shift_point_back(1,:);
    rotated_y = shift_point_back(2,:);



    figure(4);
    imagesc(distance_map);
    hold on;

    plot(rotated_x, rotated_y, 'g');
    plot(rotated_x, rotated_y, '*b');
    plot(xc, yc, 'r');
    plot(xc, yc, '*r');
    hold off;

end



if(method == 2)
    blue_point_number = length(pointx);

    [non_zero_y, non_zero_x] = find(edge_image);


    figure(3);
    imagesc(distance_map);
    hold on;
    plot(pointx, pointy, 'g');
    plot(pointx, pointy, '*b');



    for iter = 1:500
        red_point_x = zeros(1,blue_point_number);
        red_point_y = zeros(1,blue_point_number);

        for i = 1:blue_point_number
            min_blue_value = 100000000;
            min_x = 100000;
            min_y = 100000;

            for j = 1:length(non_zero_x)
                distance_squre = (pointx(i) - non_zero_x(j))^2 + (pointy(i) - non_zero_y(j))^2;
                if(distance_squre < min_blue_value)
                    min_blue_value = distance_squre;
                    min_x = non_zero_x(j);
                    min_y = non_zero_y(j);
                end        
            end

            red_point_x(i) = min_x;
            red_point_y(i) = min_y;

        end

        blue_point_x_mean = mean(pointx);
        blue_point_y_mean = mean(pointy);
        red_point_x_mean = mean(red_point_x);
        red_point_y_mean = mean(red_point_y);


        anti_mean_blue_point_x = pointx - blue_point_x_mean;
        anti_mean_blue_point_y = pointy - blue_point_y_mean;
        anti_mean_red_point_x = red_point_x - red_point_x_mean;
        anti_mean_red_point_y = red_point_y - red_point_y_mean;


        combined_blue_point = cat(1, anti_mean_blue_point_x, anti_mean_blue_point_y);
        combined_red_point = cat(1, anti_mean_red_point_x, anti_mean_red_point_y);




        H = combined_blue_point * combined_red_point';
        [U, ~, V] = svd(H);

        R = V * U';

        t = [red_point_x_mean; red_point_y_mean] - R * [blue_point_x_mean; blue_point_y_mean];

        transformed_blue_point = R * combined_blue_point + t + [blue_point_x_mean; blue_point_y_mean];

        transformed_blue_point_x = transformed_blue_point(1,:);
        transformed_blue_point_y = transformed_blue_point(2,:);
        pointx = transformed_blue_point_x;
        pointy = transformed_blue_point_y;
        plot(pointx,pointy, '*y');
        plot(pointx,pointy, 'y');
        plot(blue_point_x_mean, blue_point_y_mean,'+g', red_point_x_mean, red_point_y_mean, '+r');
    end
    plot(pointx,pointy, '*w');
    plot(pointx,pointy, 'w');
    plot(xc, yc, 'r');
    plot(xc, yc, '*r');
    hold off;

    figure(4);
    imagesc(distance_map);
    hold on;
    plot(pointx,pointy, '*w');
    plot(pointx,pointy, 'w');
    plot(xc, yc, 'r');
    plot(xc, yc, '*r');
    hold off;
    
end



function [new_xpos, new_ypos] = resample(nsample, initcontourx, initcontoury)
    len_initcontour = max(size(initcontourx));
    new_x_pos = zeros(1);
    new_y_pos = zeros(1);
    new_x_pos(2:len_initcontour) = initcontourx(2:len_initcontour) - initcontourx(1:len_initcontour-1);
    new_y_pos(2:len_initcontour) = initcontoury(2:len_initcontour) - initcontoury(1:len_initcontour-1);
    new_xx_pos_sq = new_x_pos.^2;
    new_yy_pos_sq = new_y_pos.^2;
    new_xy_pos = sqrt(new_xx_pos_sq + new_yy_pos_sq);

    cumu_sum = cumsum(new_xy_pos);

    total_length = cumu_sum(len_initcontour);
    delta = total_length/nsample;

    samples = [0:delta:total_length];

    new_xpos = interp1(cumu_sum, initcontourx, samples);
    new_ypos = interp1(cumu_sum, initcontoury, samples);
    new_xpos = new_xpos';
    new_ypos = new_ypos';
end

function process_image = create_edge(x, y, edge_image)
 process_image = edge_image;
 for i = 1:size(x)
     process_x = x(i);
     process_y = y(i);
     x_upperbound = ceil(process_x);
     x_lowerbound = floor(process_x);
     y_upperbound = ceil(process_y);
     y_lowerbound = floor(process_y);
     
     diff_x_upperbound = x_upperbound - process_x;
     diff_y_upperbound = y_upperbound - process_y;
     
     process_image(x_upperbound,y_upperbound) = process_image(x_upperbound,y_upperbound) + (1-diff_x_upperbound) * (1-diff_y_upperbound);
     process_image(x_upperbound,y_lowerbound) = process_image(x_upperbound,y_lowerbound) + (1-diff_x_upperbound) * diff_y_upperbound;
     process_image(x_lowerbound,y_upperbound) = process_image(x_lowerbound,y_upperbound) + diff_x_upperbound * (1-diff_y_upperbound);
     process_image(x_lowerbound,y_lowerbound) = process_image(x_lowerbound,y_lowerbound) + diff_x_upperbound * diff_y_upperbound; 
     
     
 end
end

function output = angle_translation(param)

pointx = evalin('base', 'pointx');
pointy = evalin('base', 'pointy');
distance_map = evalin('base', 'distance_map');
v = cat(1,pointx,pointy);


x_center = 128;
y_center = 128;

rotation_center = repmat([x_center; y_center], 1, length(pointx));
rotation_matrix = [cos(param(1)) -sin(param(1)); sin(param(1)) cos(param(1))];

rotation_center_translation = repmat([param(2); param(3)], 1, length(pointx));
shift_point = v - rotation_center;     
shift_point_applied_origion = rotation_matrix*shift_point;           
shift_point_back = shift_point_applied_origion + rotation_center + rotation_center_translation;  

rotated_x = shift_point_back(1,:);
rotated_y = shift_point_back(2,:);

output = sum(interp2(distance_map, rotated_x, rotated_y));



disp("ASSA");
plot(rotated_x, rotated_y, 'y');

end
