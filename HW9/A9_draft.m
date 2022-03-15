file = Tiff('hand2.tif','r');
img = read(file);

% I = imread('hand1.tif');
% figure(1)
% imshow(I)

figure(1)
imagesc(img) 
% colormap(gray)
hold on;
% plot(xc1, yc1, '*r');
xc1 = [];
yc1 = [];
button = 1;
count = 0;
while button == 1
    [x_coord, y_coord, button] = ginput(1);
    count = count+1;
    xc1(count) = x_coord;
    yc1(count) = y_coord;
    plot(x_coord, y_coord, '*r', x_coord, y_coord, 'r');
end
hold off;
save('h2.mat','xc1','yc1');