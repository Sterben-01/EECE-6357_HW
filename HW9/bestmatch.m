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