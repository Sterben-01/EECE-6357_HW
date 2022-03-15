function [nearx, neary, dist] = p2l(px1, px2, py1, py2, qx, qy)
    % This function returns the minimum distance between point (qx, qy) to
    % the line segment with terminals (px1, py1) and (px2, py2)
    u = ((qx - px1) * (px2 - px1) + (qy - py1) * (py2 - py1)) / ...
        ((px2 - px1)^2 + (py2 - py1)^2);
    u = max([min([u, 1]), 0]);
    nearx = px1 + u * (px2 - px1);
    neary = py1 + u * (py2 - py1);
    dist = sqrt((qx - nearx)^2 + (qy - neary)^2);
end

