function [processed_point_x, processed_point_y] = ICP (initial_point_x, initial_point_y, target_point_x, target_point_y)
    
    blue_point_number = length(initial_point_x);


    for iter = 1:500
        red_point_x = zeros(1,blue_point_number);
        red_point_y = zeros(1,blue_point_number);

        for i = 1:blue_point_number
            min_blue_value = 100000000;
            min_x = 100000;
            min_y = 100000;

            for j = 1:length(non_zero_x)
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

        transformed_blue_point = R * combined_blue_point + t + [blue_point_x_mean; blue_point_y_mean];

        transformed_blue_point_x = transformed_blue_point(1,:);
        transformed_blue_point_y = transformed_blue_point(2,:);
        initial_point_x = transformed_blue_point_x;
        initial_point_y = transformed_blue_point_y;

    end
    processed_point_x = initial_point_x;
    processed_point_y = initial_point_y;
    
end