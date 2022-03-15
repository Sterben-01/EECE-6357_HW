function [out]=snake_mainfunction(binim, forcetype, std, support, Niter, nsample, alpha, beta, gamma, extcoef, balcoef, itergvf)
    sig = std;
    figure(1);
    subplot(1,1,1),imagesc(binim);
    binim = double(binim);
    hold on;
    initcontourx = evalin('base', 'initcontourx');
    initcontoury = evalin('base', 'initcontoury');
    plot(initcontourx,initcontoury,'r');
    plot(initcontourx,initcontoury,'*r');
    hold off;
    
    if (forcetype == 1)
        [process_1_matrix_x, process_1_matrix_y] = meshgrid(-support:support);
        process_1_matrix = process_1_matrix_x.^2+process_1_matrix_y.^2;
        gaussion_value = (1/sqrt(2*pi)*sig)*exp(-(process_1_matrix./(2*sig^2)));
        process_binim = conv2(binim, gaussion_value,'same');
        binim_min = min(process_binim,[],'all');
        binim_max = max(process_binim,[],'all');
        processed_binim = (process_binim - binim_min)./(binim_max - binim_min);
        [image_forcex,image_forcey] = gradient(processed_binim);%energy method2
        forceX = image_forcex;
        forceY = image_forcey;
        figure(2);
        subplot(1,1,1),imagesc(binim), colormap('gray');
        hold on;
        quiver(forceX,forceY); %method2
        hold off;

        [added_x, added_y] = resample(nsample, initcontourx, initcontoury);

        figure(3);

        subplot(1,1,1),imagesc(binim);
        hold on;
        plot(added_x,added_y,'g');
        plot(added_x,added_y,'*g');




    end
    if (forcetype == 2)
        distforce = bwdist(binim); %energy method2
        [distforcex,distforcey] = gradient(distforce);%energy method2
        forceX = -distforcex;
        forceY = -distforcey;
        figure(2);
        subplot(1,1,1),imagesc(binim), colormap('gray');
        hold on;
        quiver(forceX,forceY); %method2
        hold off;

        [added_x, added_y] = resample(nsample, initcontourx, initcontoury);

        figure(3);

        subplot(1,1,1),imagesc(binim);
        hold on;
        plot(added_x,added_y,'g');
        plot(added_x,added_y,'*g');

    end
    if (forcetype ==3) 
        [GVFforcex, GVFforcey] = GVF(binim, 0.2, itergvf); %method 3
        forceX = GVFforcex;
        forceY = GVFforcey;
        figure(2);
        subplot(1,1,1),imagesc(binim), colormap('gray');
        hold on;
        quiver(forceX, forceY); %method 3
        hold off;
        % new_x_pos = initcontourx(2:7) - initcontourx(1:6);
        % new_y_pos = initcontoury(2:7) - initcontoury(1:6);
        % new_x_pos = [[0],new_x_pos];
        % new_y_pos = [[0],new_y_pos];
        % new_x_pos_sq = new_x_pos.*new_x_pos;
        % new_y_pos_sq = new_y_pos.*new_y_pos;
        % new_xy_pos = sqrt(new_x_pos_sq + new_y_pos_sq);
        % 
        % cumu_sum = cumsum(new_xy_pos);
        % 
        % total_length = sum(new_xy_pos);
        % delta = total_length/nsample;
        % 
        % samples = [0:delta:total_length];
        % 
        % added_x = interp1(cumu_sum, initcontourx, samples);
        % added_y = interp1(cumu_sum, initcontoury, samples);

        [added_x, added_y] = resample(nsample, initcontourx, initcontoury);

        figure(3);

        subplot(1,1,1),imagesc(binim);
        hold on;
        plot(added_x,added_y,'g');
        plot(added_x,added_y,'*g');


        % 
        % distforcex_t = distforcex;
        % distforcey_t = distforcey;
        % 
        % for i = 1:itergvf
        %     
        %     distforcex_t = 0.2*del2(binim)*4*distforcex_t-(distforcex_t-distforcex)*(distforcex_t.^2+distforcey_t.^2);
        %     distforcey_t = 0.2*del2(binim)*4*distforcey_t-(distforcey_t-distforcey)*(distforcex_t.^2+distforcey_t.^2);
        %     
        %     
        % 
        % 
        % 
        % 
        % 
        % end

    end

    alpha = alpha* ones(1,nsample+1);
    beta = beta*ones(1,nsample+1);
    % produce the five diagnal vectors
    alpham1 = [alpha(2:nsample+1) alpha(1)];
    alphap1 = [alpha(nsample+1) alpha(1:nsample)];
    betam1 = [beta(2:nsample+1) beta(1)];
    betap1 = [beta(nsample+1) beta(1:nsample)];
    a = betam1;
    b = -alpha - 2*beta - 2*betam1;
    c = alpha + alphap1 +betam1 + 4*beta + betap1;
    d = -alphap1 - 2*beta - 2*betap1;
    e = betap1;
    % generate the parameters matrix
    A = diag(a(1:nsample-1),-2) + diag(a(nsample:nsample+1),nsample-1);
    A = A + diag(b(1:nsample),-1) + diag(b(nsample+1), nsample);
    A = A + diag(c);
    A = A + diag(d(1:nsample),1) + diag(d(nsample+1),-(nsample));
    A = A + diag(e(1:nsample-1),2) + diag(e(nsample:nsample+1),-(nsample-1));
    %invAI = inv(A + gamma * diag(ones(1,nsample)));

    I = diag(ones(nsample+1,1));

    % [interp_x, interp_y] = interpForce(GVFforcex, GVFforcey, added_x, added_y);
    snakex = added_x;
    snakey = added_y;


    for i = 1:Niter

        % N is the length of snake
        size_snakex = length(snakex);
        size_snakey = length(snakey);
        snakexr = snakex;
        snakexr(2:size_snakex) = snakex(1:size_snakex - 1);
        snakexr(1) = snakex(size_snakex);
        snakexl = snakex;
        snakexl(1:size_snakex - 1) = snakex(2:size_snakex);
        snakexl(size_snakex) = snakex(1);

        snakeyr = snakey;
        snakeyr(2:size_snakey) = snakey(1:size_snakey-1);
        snakeyr(1) = snakey(size_snakey);
        snakeyl = snakey;
        snakeyl(1:size_snakey-1) = snakey(2:size_snakey);
        snakeyl(size_snakey) = snakey(1);

        balf_x_pre = (snakeyr - snakeyl) ./ sqrt((snakexr - snakexl).^2 + (snakeyr - snakeyl).^2);
        balf_y_pre = -(snakexr - snakexl) ./ sqrt((snakexr - snakexl).^2 + (snakeyr - snakeyl).^2);

        [interp_x, interp_y] = interpForce(forceX, forceY, snakex, snakey);



        snakex = (I + gamma * A) \ (snakex + extcoef*interp_x + balf_x_pre*balcoef);
        snakey = (I + gamma * A) \ (snakey + extcoef*interp_y + balf_y_pre*balcoef);
        [snakex,snakey] = resample(nsample, snakex', snakey');
        plot(snakex,snakey,'g');
        plot(snakex,snakey,'*g');
        pause(0.05);

        final_x = snakex;
        final_y = snakey;
    end

    hold off;

    figure(4);

    subplot(1,1,1);
    imagesc(binim);
    hold on;
    plot(final_x,final_y,'*g');
end


function [interp_forcex, interp_forcey] = interpForce(fieldx, fieldy, x, y)
    indx = 1:size(fieldx(1, :)');
    indy = 1:size(fieldy(1, :)');
    interp_forcex = interp2(indx, indy, fieldx, x, y);
    interp_forcey = interp2(indx, indy, fieldy, x, y);
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













