



function [calced_matrix_x, calced_matrix_y] = distance_match(shape_x, shape_y,target_x, target_y)
    %calc zhou chang
    perimeter_shape = 0;
    perimeter_target = 0;
    
    calced_matrix_x = zeros(1, length(shape_x));
    
    calced_matrix_y = zeros(1, length(target_y));
    
    
    for i = 1:length(shape_x)-1
        perimeter_shape = perimeter_shape + norm([(shape_x(i+1) - shape_x(i)), (shape_y(i+1) - shape_y(i))] ,2);
    end
    
    for i = 1:length(target_x)-1
        perimeter_target = perimeter_target + norm([(target_x(i+1) - target_x(i)), (target_y(i+1) - target_y(i))] ,2);
    end
    
    cumulate_shape_used_distance = 0;
    
    for j = 1:length(shape_x)-1
        
        temp_distance = norm([(shape_x(j+1) - shape_x(j)), (shape_y(j+1) - shape_y(j))] ,2);
        cumulate_shape_used_distance = cumulate_shape_used_distance + temp_distance;
       
        cumulate_target_used_distance = 0;
        
        
        
        for p = 2:length(target_x)
            
            temp_distance_target = norm([(target_x(p) - target_x(p-1)), (target_y(p) - target_y(p-1))] ,2);
            
            
            if((cumulate_target_used_distance + temp_distance_target) / perimeter_target) >= (cumulate_shape_used_distance / perimeter_shape)
                
                A = cumulate_shape_used_distance / perimeter_shape; 
                B = cumulate_target_used_distance / perimeter_target;
                C = (cumulate_target_used_distance + temp_distance_target) / perimeter_target;
                
                precentage_in_use = (A-B)/(C-B);
                
                delta = (cumulate_target_used_distance + temp_distance_target) / perimeter_target - cumulate_shape_used_distance / perimeter_shape;
                portion = 1-delta / (temp_distance_target / perimeter_target);
                calced_matrix_x(j) = target_x(p-1) + portion*(target_x(p) - target_x(p-1));
                calced_matrix_y(j) = target_y(p-1) + portion*(target_y(p) - target_y(p-1));
                
                break;
                
            end
            
            cumulate_target_used_distance = cumulate_target_used_distance + temp_distance_target;
            
        end

    end
    
    calced_matrix_x(1) = target_x(1);
    calced_matrix_y(1) = target_y(1);
    calced_matrix_x(length(shape_x)) = target_x(length(target_x));
    calced_matrix_y(length(shape_x)) = target_y(length(target_y));
    
end


% shape 69x2
% target n*2
%calced_matrix 69x2

function [x, y] = bestmatch(shape, target)
    
    n_shape = length(shape);
    circ_shape = 0;
    for i = 2:n_shape
        circ_shape = circ_shape + norm([shape(i, 1) - shape(i - 1, 1), ...
                                        shape(i, 2) - shape(i - 1, 2)], 2);
    end
    
    n_target = length(target);
    circ_target = 0;
    for i = 2:n_target
        circ_target = circ_target + norm([target(i, 1) - target(i - 1, 1), ...
                                          target(i, 2) - target(i - 1, 2)], 2);
    end
    
    x = zeros(1, n_shape);
    y = zeros(1, n_shape);
    length_shape = 0;
    for i = 1:n_shape - 1
        length_target = 0;
        if i > 1
            length_shape = length_shape + norm([shape(i, 1) - shape(i - 1, 1), ...
                                                shape(i, 2) - shape(i - 1, 2)], 2);
        end
        for j = 2:n_target
            len = norm([target(j, 1) - target(j - 1, 1), ...
                        target(j, 2) - target(j - 1, 2)], 2);
            if (length_target + len) / circ_target <  length_shape / circ_shape
                length_target = length_target + len;
            else
                delta = (length_target + len) / circ_target - length_shape / circ_shape;
                portion = 1 - delta / (len / circ_target);
                x(i) = target(j - 1, 1) + portion * (target(j, 1) - target(j - 1, 1));
                y(i) = target(j - 1, 2) + portion * (target(j, 2) - target(j - 1, 2));
                break
            end
        end
    end
    
    x(1) = target(1, 1);
    y(1) = target(1, 2);
    x(n_shape) = target(n_target, 1);
    y(n_shape) = target(n_target, 2);
end






