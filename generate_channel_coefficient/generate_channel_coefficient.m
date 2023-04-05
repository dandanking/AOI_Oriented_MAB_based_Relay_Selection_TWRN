%%% 随机生成信道系数(20*4):
%%% 可划分为10组。每组2*4


clc
clear all

disp('-----------------------')
disp('生成10组随机信道系数')
disp('-----------------------')

% Random_channel_coefficient = zeros(10,8);
random_channel_coefficient = 1/sqrt(2)*(randn(20,4) + j*randn(20,4));

save('channel_coefficient_(20_4)_1.mat');
