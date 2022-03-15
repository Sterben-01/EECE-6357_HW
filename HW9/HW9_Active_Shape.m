%% Phase 1: hand registration and getting mean shape

X = zeros(10, 64);
Y = zeros(10, 64);

figure(1);
subplot(2, 2, 1);
for i = 1:10
    X(i, :) = load(sprintf('h%d.mat', i)).x(:);
    Y(i, :) = load(sprintf('h%d.mat', i)).y(:);
    plot(X(i, :), Y(i, :), 'b*', X(i, :), Y(i, :), 'b-');
    hold on;
end
hold off;
title('Initial shapes');

de = zeros(10, 64, 2);
for i = 1:10
    de(i, :, 1) = X(i, :);
    de(i, :, 2) = Y(i, :);
end

for k = 1:10
    mean_shape = squeeze(de(1, :, :));
    % Register hands 2~10 to hand 1
    for j = 2:10
        [dx, dy, ~, ~] = ICP(de(j, :, 1), de(j, :, 2), de(1, :, 1), de(1, :, 2));
        de(j, :, 1) = dx(:);
        de(j, :, 2) = dy(:);
    end
    new_mean_shape = squeeze(mean(de, 1));
    if norm(mean_shape - new_mean_shape) < 0.1
        break
    end
    de(1, :, :) = new_mean_shape;
end

subplot(2, 2, 2);
for i = 1:10
    plot(squeeze(de(i, :, 1)), squeeze(de(i, :, 2)), 'b*', ...
         squeeze(de(i, :, 1)), squeeze(de(i, :, 2)), 'b-');
    hold on;
end
hold off;
xlim([0 450]);
ylim([0 450]);
title('Registered shapes');

subplot(2, 2, 3);
mean_shape = squeeze(mean(de, 1));
plot(mean_shape(:, 1)', mean_shape(:, 2)', 'r-');
xlim([0 450]);
ylim([0 450]);
title('Mean shape');

%% Phase 2: Eigenshape
features = reshape(de - repmat(reshape(mean_shape, 1, 64, 2), 10, 1, 1), 10, 128);
cov_mat = cov(features);
[evec, ~, explained] = pcacov(cov_mat);

explained = explained ./ sum(explained);
useful_index = find(cumsum(explained) <= 0.95);
used_index = max(useful_index);
eigen_weight = evec(:, 1:used_index);

figure(2);
for i = 1:4
    subplot(2, 2, i);
    for factor = -100:10:100
        shape_disp = mean_shape + reshape(factor * evec(:, i), 64, 2);
        plot(shape_disp(:, 1), shape_disp(:, 2), 'b-');
        xlim([0 400]);
        ylim([0 450]);
        title(sprintf('eigen shape %d', i));
        pause(0.1);
    end
end


%% Phase 3: Fitting to new images

% Read in the data
while true
    try
        in_name = input('The name of your test file is: ', 's');
        in_name = ProcessFilename(in_name);
        try
            xy = load(in_name).xy;
        catch
            xy = load(in_name).xy2;
        end
        break
    catch
        fprintf('File does not exist, please input again.\n');
    end
end

% Initialize b = 0
b = zeros(used_index, 1);

[n_targets, ~] = size(xy);
[n_fid, ~] = size(mean_shape);
figure(3);
subplot(2, 2, 1);
plot(xy(:, 1), xy(:, 2), 'g*', xy(:, 1), xy(:, 2), 'g-');
hold on;
plot(mean_shape(:, 1), mean_shape(:, 2), 'r-');
hold off;
title('Initial state');

subplot(2, 2, 2);
% Iteration
for i = 1:20
    % Generate current shape model
    shape = mean_shape + reshape(eigen_weight * b, n_fid, 2);
    x_shape = shape(:, 1)';
    y_shape = shape(:, 2)';
    x = xy(:, 1)';
    y = xy(:, 2)';
    [x_disp, y_disp, T, R] = ICP(x_shape, y_shape, x, y);
    
    plot(x, y, 'g*-');
    hold on;
    plot(x_disp, y_disp, 'b*-');
    hold off;
    xlim([0 450]);
    ylim([0 450]);
    title(sprintf('Iter %d', i));

    [nx, ny] = bestmatch(shape, xy);
    
    xy_nearest = cat(1, nx, ny)';
    xy_inverted = (R \ (xy_nearest - repmat(T, 1, n_fid)')')';  % * inv(R);
    x_inverted = xy_inverted(:, 1)';
    y_inverted = xy_inverted(:, 2)';
    
    b = eigen_weight' * reshape((cat(1, x_inverted, y_inverted)' - mean_shape), 2 * n_fid, 1);

    pause(0.5);
end
hold off;
xlim([0 450]);
ylim([0 450]);

shape = mean_shape + reshape(eigen_weight * b, n_fid, 2);
x_shape = shape(:, 1)';
y_shape = shape(:, 2)';
[x_shape, y_shape, T, R] = ICP(x_shape, y_shape, x, y);

subplot(2, 2, 3);

plot(x, y, 'g*', x, y, 'g-');
hold on;
plot(x_shape, y_shape, 'b*', x_shape, y_shape, 'b-');
hold off;
xlim([0 450]);
ylim([0 450]);
title('Final results');

