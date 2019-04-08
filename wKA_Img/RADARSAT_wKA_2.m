% RADARSAT-1的数据处理之wKA成像2
% wKA的近似实现

clc;clear all;close all;

%% 提取原始数据
run specify_parameters.m
run extract_data.m
clear all;

%% 参数设置
% 雷达参数
C = 2.9979e8;                                                               % 光速
f0 = 5.3e9;                                                                 % 中心频率
lambda = C/f0;                                                              % 波长
% 飞行平台
V = 7062;                                                                   % 平台速度
% 距离向
Kr = -0.72135e12;                                                           % 调频率
Tr = 41.75e-6;                                                              % 脉冲持续时间
Fr = 32.317e6;                                                              % 距离向采样率
dtau = 1/Fr;                                                                % 距离向采样时间间隔
% 方位向
Ka = 1733;                                                                  % 方位向调频率
Fa = 1256.98;                                                               % 方位向采样率
fc = -6900;                                                                 % 多普勒中心频率
deta = 1/Fa;                                                                % 距离向采样时间间隔
% 其他参数
t0 = 6.5956e-3;                                                             % 获取数据时间
R0 = t0*C/2;                                                                % 最短斜距

%% 获取原始数据
oldFolder = cd('.\scene01');
load CDdata1.mat
cd(oldFolder);
data=double(data);
[Na,Nr] = size(data);

%% 确定参考距离
Rref = (t0+Nr/2*dtau)*C/2;                                                  % 参考距离选在成像区域中心

%% 数据末尾补零
Za = 800;                                                                   % 方位向补零数
Zr = ceil(Tr/dtau);                                                         % 距离向补零数
data = cat(2,data,zeros(Na,Zr));                                            % 距离向补零
data = cat(1,zeros(Za,Nr+Zr),data);                                         % 方位向补零
Na = Na+Za;
Nr = Nr+Zr;

figure,imagesc(abs(data));axis image;set(gcf,'Color','w');
title('时域：补零后的原始信号幅度');
xlabel('距离向（采样点）');ylabel('方位向（采样点）');

%% 时间轴、频率轴设置
tau = t0+(0:Nr-1)*dtau;                                                     % 距离向时间轴
ftau = ((0:Nr-1)-Nr/2)/Nr*Fr;                                               % 距离向频率轴
eta = ((0:Na-1)-Na/2)*deta;                                                 % 方位向时间轴
feta = fc+((0:Na-1)-Na/2)/Na*Fa;                                            % 方位向频率轴

%% 对齐到方位向多普勒中心
S0 = data.*exp(-1i*2*pi*fc*(eta'*ones(1,Nr)));                              % 搬移至多普勒频率中心
clear data

%% 参考函数相乘
Sff = fftshift(fft2(fftshift(S0)));                                         % 二维FFT
clear S0
WindowR = ones(Na,1)*kaiser(Nr,2.5)';                                       % 距离窗
WindowA = kaiser(Na,2.5)*ones(1,Nr);                                        % 方位窗
Sff = Sff.*WindowR.*WindowA;                                                % 加窗
clear WindowR WindowA
Hrfm = exp(1i*pi*(4*Rref/C*sqrt((f0+ones(Na,1)*ftau).^2-(C*(feta'*ones(1,Nr))/V/2).^2)+(ones(Na,1)*ftau).^2/Kr)); % 构造参考函数，式8.3
S1 = Sff.*Hrfm;                                                             % 参考函数相乘（一致压缩），式8.1
clear Sff Hrfm

figure,imagesc(abs(S1));axis image;set(gcf,'Color','w');
title('二维频域：参考函数相乘后的信号幅度');
xlabel('距离向（采样点）');ylabel('方位向（采样点）');

%% 补余方位匹配滤波器相乘
Srd = ifftshift(ifft(ifftshift(S1,2),[],2),2);                              % 距离向IFFT
clear S1
D = sqrt(1-C^2*feta.^2/(4*V^2*f0^2));                                       % 距离徙动因子D(feta,V)，式7.17
Happ = exp(1i*4*pi*((ones(Na,1)*tau)*C/2-Rref)*f0.*(D'*ones(1,Nr))/C);      % 补余方位匹配滤波器，式8.35指数项相消，R0?
clear D
S2 = Srd.*Happ;                                                             % 滤波相乘
clear Srd Happ

figure,imagesc(abs(S2));axis image;set(gcf,'Color','w');
title('RD域：补余方位匹配滤波器相乘后的信号幅度');
xlabel('距离向（采样点）');ylabel('方位向（采样点）');

%% 变换到图像域
Stt = ifftshift(ifft(ifftshift(S2,1),[],1),1);                              % 方位向IFFT
clear S2

%% 调整图像
Img = flipud(abs(Stt));                                                     % 上下翻转图像
Ntau = round((Tr/dtau/2-1)/2);                                              % 计算距离向弃置区长度
Img = Img(1:Na-Za,Ntau+1:Ntau+Nr-Zr);                                       % 裁剪图像有效区域
Img = Img/max(max(Img));
Img = 20*log10(Img+eps);
Img(Img<-50) = -50;
Img = uint8((Img+50)/50*255);

figure,imagesc(Img);axis image;set(gcf,'Color','w');
figure,imshow(Img);
imwrite(Img,'SAR_wKA_2.bmp','bmp');
