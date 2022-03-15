file = input('please input image name: ', 's');
load(file);
process = input('please input process mode: ');
hist = imhist(testima); %create histogram
L = 0:255; %gray color range
len = length(hist);
prob = hist./sum(hist); %normalization
param = input('please input parameter: ');
%param = [0.3, 10, 60, 0.4, 10, 100, 0.3, 10, 150];
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

reshape_img = reshape(double(testima), 1, 65536);
if (process == 0)
    figure(1); %create figure
    subplot(2,2,1), plot(L, hist); % create 2*2 skeleton No.1 graph
        
        
    %normal probability density function
    class_1 = prob_1*(1/sqrt(2*pi*sig_1^2))*exp(-(mu_1-L).^2/(sig_1^2*2));
    class_2 = prob_2*(1/sqrt(2*pi*sig_2^2))*exp(-(mu_2-L).^2/(sig_2^2*2));
    class_3 = prob_3*(1/sqrt(2*pi*sig_3^2))*exp(-(mu_3-L).^2/(sig_3^2*2));
    norm_dis = class_1+class_2+class_3;
    norm_dis = norm_dis/sum(norm_dis);

    subplot(2,2,2), plot(norm_dis); % Initial distribution

    options = optimset('MaxIter', 1000, 'TolFun', 1e-2, 'TolX', 1e-2);
    optim = fminsearch(@error, param, options);

    %subplot(2,2,3);
    %plot(L, prob, 'b', L, norm_dis, 'r');
    clac_1 = optim(1)*(1/sqrt(2*pi*optim(2)^2))*exp(-(optim(3)-L).^2/(optim(2)^2*2));
    clac_2 = optim(4)*(1/sqrt(2*pi*optim(5)^2))*exp(-(optim(6)-L).^2/(optim(5)^2*2));
    clac_3 = optim(7)*(1/sqrt(2*pi*optim(8)^2))*exp(-(optim(9)-L).^2/(optim(8)^2*2));
    after_calc = clac_1+clac_2+clac_3;
    after_calc = after_calc/sum(after_calc);
    subplot(2,2,4);
    plot(L, after_calc, 'b', L, prob, 'r');

elseif (process == 1)
    prob = prob';
    figure(1); %create figure
    subplot(2,2,1), plot(L, prob); % create 2*2 skeleton No.1 graph
    fxhat = prob_1*(1/sqrt(2*pi*sig_1^2))*exp(-(mu_1-L).^2/((sig_1^2)*2))+prob_2*(1/sqrt(2*pi*sig_2^2))*exp(-(mu_2-L).^2/((sig_2^2)*2))+prob_3*(1/sqrt(2*pi*sig_3^2))*exp(-(mu_3-L).^2/((sig_3^2)*2));
    subplot(2,2,2), plot(fxhat);
    for LL = 0:200
    subplot(2,2,3);
        
    %plot(L, hist);
    %hold on;
        
    pause(0.1);
        
    %plot(fxhat, 'b');
    %hold off;
        
    plot(L, prob, 'b', L, fxhat, 'r');
        
    class_1 = (1/sqrt(2*pi*sig_1^2))*exp(-(mu_1-reshape_img).^2/((sig_1^2)*2));
    class_2 = (1/sqrt(2*pi*sig_2^2))*exp(-(mu_2-reshape_img).^2/((sig_2^2)*2));

        
    %disp(class_1);
    %disp(class_2);
        
    %single distripution occurance probability
        
        
    w1 = class_1*prob_1./(class_1*prob_1+class_2*prob_2);
    w2 = class_2*prob_2./(class_1*prob_1+class_2*prob_2);
    
        
    %disp(w1);
    %disp(w2);

        
    %recalculate distripution
        
    %prob_1 = prob*w1';
    prob_1 = (1/65536)*sum(w1);
    prob_2 = (1/65536)*sum(w2);

        
    %disp(prob_1);
    %disp(prob_2);
        
    %mean
        
    mu_1 = (1/(65536*prob_1))*sum(w1.*reshape_img);
    mu_2 = (1/(65536*prob_2))*sum(w2.*reshape_img);
    
        
    %disp(mu_1);
    %disp(mu_2);
        
    % variance
        
    sig_1 = (1/(65536*prob_1)*sum((reshape_img-mu_1).^2.*w1));
    sig_2 = (1/(65536*prob_2)*sum((reshape_img-mu_2).^2.*w2));

        
    %disp(sig_1);
    %disp(sig_2);
      
    %standard devation
        
    sig_1 = sqrt(sig_1);
    sig_2 = sqrt(sig_2);
        
    %disp(sig_1);
    %disp(sig_2);

        
        
    %correction image gaussion distribution.
        
    fxhat = prob_1*(1/sqrt(2*pi*sig_1^2))*exp(-(mu_1-L).^2/((sig_1^2)*2))+prob_2*(1/sqrt(2*pi*sig_2^2))*exp(-(mu_2-L).^2/((sig_2^2)*2));
        
        
    %subplot(2,2,3);
        
    %plot(L, hist);
    %hold on;
        
    %pause(0.01);
        
    %plot(fxhat, 'b');
    %hold off;
        
    %plot(L, prob, 'b', L, fxhat, 'r');
        
        
end
    subplot(2,2,4);
    plot(L, prob, 'k', L, fxhat, 'g+');
end






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
function err = error(param)
hist = evalin('base', 'hist');
L = 0:255;
class_A = param(1)*(1/sqrt(2*pi*param(2)^2))*exp(-(param(3)-L).^2/(param(2)^2*2));
class_B = param(4)*(1/sqrt(2*pi*param(5)^2))*exp(-(param(6)-L).^2/(param(5)^2*2));
class_C = param(7)*(1/sqrt(2*pi*param(8)^2))*exp(-(param(9)-L).^2/(param(8)^2*2));
norm_dis_func = class_A+class_B+class_C;
norm_dis_func = norm_dis_func/sum(norm_dis_func);
prob = hist./sum(hist);
prob = prob';
err = sum((prob-norm_dis_func).^2);
subplot(2,2,3);
pause(0.001);
cla;
plot(L, prob, 'b', L, norm_dis_func, 'r');
end








