function [px_iter, py_iter, T, R] = ICP(px, py, oox, ooy)
    num_fid = length(px);
    px_iter = px;
    py_iter = py;
    
    T = [0; 0];
    R = [1, 0; 0, 1];
    
    for iter = 1:100
        % Find nearest points for each fiducial point
        nx = zeros(1, num_fid);
        ny = zeros(1, num_fid);
        for pid = 1:num_fid
            min_dist = 999999;
            min_x = 0;
            min_y = 0;
            for j = 1:length(oox)-1
                [nearx, neary, cur_dist] = p2l(oox(j), oox(j + 1), ooy(j), ...
                    ooy(j + 1), px_iter(pid), py_iter(pid));
                if cur_dist < min_dist
                    min_dist = cur_dist;
                    min_x = nearx;
                    min_y = neary;
                end
            end
            nx(pid) = min_x;
            ny(pid) = min_y;
        end

        % Demean fiducial and target
        mean_px = mean(px_iter);
        mean_py = mean(py_iter);
        mean_nx = mean(nx);
        mean_ny = mean(ny);
        demean_px = px_iter - mean_px;
        demean_py = py_iter - mean_py;
        demean_nx = nx - mean_nx;
        demean_ny = ny - mean_ny;

        % Derive rotation matrix and translation vector
        H = cat(1, demean_px, demean_py) * cat(1, demean_nx, demean_ny)';
        [U, ~, V] = svd(H);
        rotation_mat = V * U';
        trans = [mean_nx; mean_ny] - rotation_mat * [mean_px; mean_py];

        % Obtain the new fiducial points
        transformed_coord = trans + rotation_mat * cat(1, px_iter, py_iter);
        px_iter_new = transformed_coord(1, :);
        py_iter_new = transformed_coord(2, :);

        T = T + trans;
        R = R * rotation_mat;
        
        if norm(px_iter - px_iter_new) < 1e-3 && ...
           norm(py_iter - py_iter_new) < 1e-3
            break
        end
        px_iter = px_iter_new;
        py_iter = py_iter_new;
    end
end

