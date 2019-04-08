% RADARSAT-1的数据处理之wKA成像1
% wKA的精确实现

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

%% Stolt插值
% 频率轴映射
fnew_m = sqrt((f0+ones(Na,1)*ftau).^2-(C*(feta'*ones(1,Nr))/V/2).^2)-f0;    % 新的频率矩阵，式8.5
fmin = min(min(fnew_m));                                                    % 频率最小值
fmax = max(max(fnew_m));                                                    % 频率最大值
fnew_m = ones(Na,1)*linspace(fmin,fmax,Nr);                                 % 更新新的频率矩阵
fold_m = sqrt((f0+fnew_m).^2+(C*(feta'*ones(1,Nr))/V/2).^2)-f0;             % 由新的频率，计算旧的频率
clear fnew_m
% 插值
Sstolt = zeros(Na,Nr);                                                      % 新建矩阵，存放插值后的数据
h = waitbar(0,'Stolt插值');
tic
for ii = 1:Na
    for jj = 1:Nr
        Delta = (fold_m(ii,jj)-ftau(jj))/(Fr/Nr);                           % 计算相对频率差，并量化
        IntNum = ceil(Delta);                                               % 频率差向上取整
        kk = jj+IntNum;                                                     % 在原矩阵中的索引
        if(5<=kk && kk<=Nr-3)                                               % 边界限定
            DecNum = IntNum-Delta;                                          % 频率差的小数部分
            SincVal = sinc(DecNum-4:1:DecNum+3)';                           % 生成插值核
            Sstolt(ii,jj) = S1(ii,kk-4:kk+3)*SincVal;                       % 插值
        end
    end
    waitbar(ii/Na);
end
toc
close(h);
clear S1 fold_m

%% 变换到图像域
Stt = ifftshift(ifft2(ifftshift(Sstolt)));                                  % 二维IFFT
% clear Sstolt

%% 调整图像
Img = flipud(abs(Stt));                                                     % 上下翻转图像
Img = Img(1:Na-Za,1:Nr-Zr);                                                 % 裁剪图像有效区域
Img = Img/max(max(Img));
Img = 20*log10(Img+eps);
Img(Img<-50) = -50;
Img = uint8((Img+50)/50*255);

figure,imagesc(Img);axis image;set(gcf,'Color','w');
figure,imshow(Img);
imwrite(Img,'SAR_wKA_1.bmp','bmp');
