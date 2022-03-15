close all;
clc;
file = input('please input image name: ', 's');
load(file);
hist = imhist(ima_att); %create histogram
L = 0:255; %gray color range
len = length(hist);
prob = hist./sum(hist); %normalization
figure(1); %create figure
set(gcf,'Position',[500,500,1900,500]);
subplot(1,3,1), surf(L, L,ima_att,'EdgeColor','flat'); % original image 3D

subplot(1,3,2), plot(L, prob); %original distribution

subplot(1,3,3), imagesc(L, L,ima_att); % original image 2D


prob_1 = param(3);
sig_1 = param(2);
mu_1 = param(1);        
prob_2 = param(6);
sig_2 = param(5);
mu_2 = param(4);


dis_scale = ima_att(1:8:256, 1:8:256);

reshape_img = double(reshape(dis_scale, 1, 1024));

matrix_phi = zeros(1024,6);

%a,b,c,d,e,f,n = 0; %ax^2 + by^2 + cxy + dx + ey + f
a = 0;
b = 0;
c = 0;
d = 0;
e = 0;
f = 0;
n = 1;

param_vec = [a,b,c,d,e,f]; % let param become a vector


%create phi matrix
for x = 0:8:255
    for y = 0:8:255
        
        phi_1 = [x^2, y^2, x*y, x, y, 1];
        
        matrix_phi(n,:) = phi_1;
        
        n = n+1;
        
        
    end
end



for LL = 0:200
    
    
    bias_func = sum(matrix_phi.*param_vec,2); %bias area
    bias_func = bias_func';    
    class_1 = (1/(sqrt(2*pi*sig_1^2)))*exp(-(reshape_img-mu_1-bias_func).^2/((sig_1^2)*2)); %bias PDF
    class_2 = (1/(sqrt(2*pi*sig_2^2)))*exp(-(reshape_img-mu_2-bias_func).^2/((sig_2^2)*2));

        
    %disp(class_1);
    %disp(class_2);
        
    %occurance probability
        
        
    w1 = class_1*prob_1./(class_1*prob_1+class_2*prob_2);
    w2 = class_2*prob_2./(class_1*prob_1+class_2*prob_2);
    
        
    %disp(w1);
    %disp(w2);

        
    %ratio 
        
    %prob_1 = prob*w1';
    prob_1 = sum(w1)*(1/1024);
    prob_2 = sum(w2)*(1/1024);

        
    %disp(prob_1);
    %disp(prob_2);
        
    %mean with bias
        
    mu_1 = sum((reshape_img-bias_func).*w1)*(1/(1024*prob_1));
    mu_2 = sum((reshape_img-bias_func).*w2)*(1/(1024*prob_2));
    
        
    %disp(mu_1);
    %disp(mu_2);
        
    %standard devation
        
    sig_1 = (sum(((reshape_img-mu_1-bias_func).^2).*w1))*(1/(1024*prob_1));
    sig_2 = (sum(((reshape_img-mu_2-bias_func).^2).*w2))*(1/(1024*prob_2));  
    
    %disp(sig_1);
    %disp(sig_2);
      
    
        
    sig_1 = sqrt(sig_1);
    sig_2 = sqrt(sig_2);
        
    %disp(sig_1);
    %disp(sig_2);

    alpha_1 = w1./(sig_1^2); %alpha
    alpha_2 = w2./(sig_2^2);
        
    %disp(alpha_1);
    %disp(alpha_2);
    
    noise_surface = reshape_img - ((alpha_1.*mu_1)+(alpha_2.*mu_2))./(alpha_1+alpha_2);
    noise_surface = noise_surface';
    omg_diag = diag(alpha_1+alpha_2); % create diagonal matrix 
    
    param_vec = (matrix_phi'*omg_diag)*matrix_phi;
    
    %param_vec = inv(param_vec); inverse use this way will slow and less precisely.
   
    param_vec = param_vec\((matrix_phi'*omg_diag)*noise_surface); % final parameter
    param_vec = param_vec';
    
    gauss_esti = prob_1*(1/sqrt(2*pi*sig_1))*exp(-(mu_1-L).^2/((sig_1^2)*2))+prob_2*(1/sqrt(2*pi*sig_2))*exp(-(mu_2-L).^2/((sig_2^2)*2));
    gauss_esti = gauss_esti/sum(gauss_esti);
    
        
end

surface_estimate = zeros(256,256); %create the estimate surfase as a 256 matrix;

for i = 0:255
    for j = 0:255
        phi_2 = [i^2, j^2, i*j, i, j, 1];
        surface_estimate(i+1, j+1) = sum(phi_2.*param_vec);
    
    end
end

corrected = uint8(double(ima_att)-surface_estimate');


figure(2);
subplot(3,2,1), plot(L,ima_att(:,128)); % original 0-128

subplot(3,2,2), plot(L,corrected(:,128)); % correct 0-128

subplot(3,2,3), plot(L,ima_att(128,:)); % original 128-

subplot(3,2,4), plot(L,corrected(128,:)); % correct 128-

subplot(3,2,5), plot(L,surface_estimate(:,128)); % correction surface 0-128

subplot(3,2,6), plot(L,surface_estimate(128,:)); % correction surface 128-

figure(3)
subplot(2,2,1), plot(L, prob, L, gauss_esti);

subplot(2,2,2),surf(L, L,surface_estimate,'EdgeColor','flat'); % surface_estimate 3D

subplot(2,2,3), imagesc(L, L,ima_att); % original image 2D

subplot(2,2,4), imagesc(L, L,corrected); % final image 2D

