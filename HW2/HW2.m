file = input('please input image name: ', 's');
load(file);
hist = imhist(testima); %create histogram
L = 0:255; 
len = length(hist);
prob = hist./sum(hist); %normalization
prob = prob';
figure(1); %create figure
subplot(2,2,1), plot(L, hist); % create 2*2 skeleton No.1 graph
param = [0.3, 10, 60, 0.4, 10, 140, 0.3, 10, 220];
%input section below
prob_1 = param(1);
sig_1 = param(2);
mu_1 = param(3);
prob_2 = param(4);
sig_2 = param(5);
mu_2 = param(6);
prob_3 = param(7);
sig_3 = param(8);
mu_3 = param(9);
%normal distribution function
class_1 = prob_1*(1/sqrt(2*pi*sig_1^2))*exp(-(mu_1-L).^2/(sig_1^2*2));
class_2 = prob_2*(1/sqrt(2*pi*sig_2^2))*exp(-(mu_2-L).^2/(sig_2^2*2));
class_3 = prob_3*(1/sqrt(2*pi*sig_3^2))*exp(-(mu_3-L).^2/(sig_3^2*2));
norm_dis = class_1+class_2+class_3;
norm_dis = norm_dis/sum(norm_dis);

subplot(2,2,2), plot(norm_dis); % Initial distribution

options = optimset('MaxIter', 2000, 'TolX', 1e-2);
optim = fminsearch(@error, param, options);
clac_1 = optim(1)*(1/sqrt(2*pi*optim(2)^2))*exp(-(optim(3)-L).^2/(optim(2)^2*2));
clac_2 = optim(4)*(1/sqrt(2*pi*optim(5)^2))*exp(-(optim(6)-L).^2/(optim(5)^2*2));
clac_3 = optim(7)*(1/sqrt(2*pi*optim(8)^2))*exp(-(optim(9)-L).^2/(optim(8)^2*2));
after_calc = clac_1+clac_2+clac_3;
after_calc = after_calc/sum(after_calc);
subplot(2,2,4);
plot(L, after_calc, 'b', L, prob, 'r');


%function err = error(~)
    %hist = evalin('base', 'hist');
    %L = 0:255;
    %norm_dis = evalin('base', 'norm_dis');
    %prob = evalin('base', 'prob');
    %prob = prob';
    %err = sum((prob-norm_dis).^2);
    %subplot(2,2,3);
    %pause(0.001);
    %cla;
    %plot(L, prob, 'b', L, norm_dis, 'r');
%end







