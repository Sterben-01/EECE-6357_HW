close all;
clc;
file = input('please input image name: ', 's');
load(file);

original_img = load(file).mitest;
original_img_rotate = load(file).mitestrot;
image_size_x = 256;
image_size_y = 256;
compressed_img = zeros(image_size_x:image_size_y);
compressed_img_rotate = zeros(image_size_x:image_size_y);
figure(1);
subplot(2,2,1);
imagesc(original_img);
colormap(gray);
subplot(2,2,2);
imagesc(original_img_rotate);
colormap(gray);
% calc comptessed number
compressor_num = image_size_x / nbins;






for i = 1:image_size_x
    for j = 1:image_size_y
        if(original_img(i,j) < maxint && original_img(i,j) > minint)
            compressed_img(i,j) = floor(original_img(i,j)/compressor_num); 
        else
            compressed_img(i,j) = -1;

        end
        
        if(original_img_rotate(i,j) < maxint && original_img_rotate(i,j) > minint)
            compressed_img_rotate(i,j) = floor(original_img_rotate(i,j)/compressor_num);
        else
            compressed_img_rotate(i,j) = -1;
        end

    end
end

for i = 1:image_size_x
    for j = 1:image_size_y
        compressed_img_rotate(i,j) = compressed_img_rotate(i,j) + 1;
        compressed_img(i,j) = compressed_img(i,j) + 1;
    end
end






param = 0;
options = optimset('MaxIter', 10000, 'TolFun', 1e-11, 'TolX', 1e-11);
MI = fminsearch(@calculateMI,param,options);


finalimg = imrotate(compressed_img_rotate, MI*100000,'nearest','crop');


joint_hist = zeros(nbins+1:nbins+1);


for i = 1:image_size_x
    for j = 1:image_size_y
        if(compressed_img(i,j) ~= 0 && finalimg(i,j) ~= 0)
            compressed_img_gray = compressed_img(i,j) + 1;
            compressed_img_rotate_gray = finalimg(i,j) + 1;
            joint_hist(compressed_img_gray,compressed_img_rotate_gray) = joint_hist(compressed_img_gray,compressed_img_rotate_gray) + 1; 
        end
    end
end


sum_joint_hist = sum(sum(joint_hist));

pdf_joint = joint_hist./sum_joint_hist;
subplot(2,2,4);
imagesc(pdf_joint);

%rotate b
%recalculate joint hist
%use joint hist to calc HA HB HAB
%use HA HB HAB calc MI

function calcMI = calculateMI(rotate_angle)
    image_size_x = evalin('base', 'image_size_x');
    image_size_y = evalin('base', 'image_size_y');
    compressed_img= evalin('base', 'compressed_img');
    compressed_img_rotate = evalin('base', 'compressed_img_rotate');
    nbins = evalin('base', 'nbins');
    compressed_img_rotate = imrotate(compressed_img_rotate, rotate_angle*100000,'nearest','crop');
    
    

    joint_hist = zeros(nbins+1:nbins+1);



    for i = 1:image_size_x
        for j = 1:image_size_y
            if(compressed_img(i,j) ~= 0 && compressed_img_rotate(i,j) ~= 0)
                compressed_img_gray = compressed_img(i,j) + 1;
                compressed_img_rotate_gray = compressed_img_rotate(i,j) + 1;
                joint_hist(compressed_img_gray,compressed_img_rotate_gray) = joint_hist(compressed_img_gray,compressed_img_rotate_gray) + 1; 
            end
        end
    end


    sum_joint_hist = sum(sum(joint_hist));

    pdf_joint = joint_hist./sum_joint_hist;

%   H = -sum_joint_hist .* log(pdf_joint+1e-12);
%     HA = sum(sum(pdf_joint,1) .* sum(log(pdf_joint+1e-12),1));
%     HB = sum(sum(pdf_joint,2) .* sum(log(pdf_joint+1e-12),2));
    HA = sum(sum(pdf_joint,1) .* log(sum(pdf_joint+1e-12,1)));
    HB = sum(sum(pdf_joint,2) .* log(sum(pdf_joint+1e-12,2)));
    HAB = sum(sum(pdf_joint.*log(pdf_joint+1e-12)));
    MI = HA + HB - HAB;
    calcMI = MI;
    subplot(2,2,3);
    imagesc(compressed_img_rotate);
    pause(0.1);
    subplot(2,2,4);
    imagesc(pdf_joint);
end















