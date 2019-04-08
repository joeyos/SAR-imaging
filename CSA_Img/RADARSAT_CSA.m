% RADARSAT-1的数据处理之CSA成像

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

%% 参考数据设置
Rref = (t0+Nr/2*dtau)*C/2;                                                  % 参考距离
Vref = V;                                                                   % 参考速度
fref = fc;                                                                  % 参考频率

%% 数据末尾补零
Za = 800;                                                                   % 方位向补零数Tsar*Fa
Zr = ceil(Tr/dtau);                                                         % 距离向补零数Tr*Fr
data = cat(2,data,zeros(Na,Zr));                                            % 距离向补零
data = cat(1,zeros(Za,Nr+Zr),data);                                         % 方位向补零
Na = Na+Za;%2336
Nr = Nr+Zr;%3398

figure,imagesc(abs(data));axis image;set(gcf,'Color','w');
title('时域：补零后的原始信号幅度');
xlabel('距离向（采样点）');ylabel('方位向（采样点）');

%% 时间轴、频率轴设置
tau = t0+(0:Nr-1)*dtau;                                                     % 距离向时间轴
ftau = ((0:Nr-1)-Nr/2)/Nr*Fr;                                               % 距离向频率轴
eta = ((0:Na-1)-Na/2)*deta;                                                 % 方位向时间轴
feta = fc+((0:Na-1)-Na/2)/Na*Fa;                                            % 方位向频率轴
%% 中间变量设置
D_feta_V = sqrt(1-C^2*feta.^2/(4*V^2*f0^2));                                % D(feta,V)，式7.17
D_feta_Vref = sqrt(1-C^2*feta.^2/(4*Vref^2*f0^2));                          % D(feta,Vref)，式7.17
D_fref_V = sqrt(1-C^2*fref^2/(4*V^2*f0^2));                                 % D(fref,V)，式7.17
D_fref_Vref = sqrt(1-C^2*fref^2/(4*Vref^2*f0^2));                           % D(fref,Vref)，式7.17

Ksrc = (2*V^2*f0^3*D_feta_V.^3)./(C*R0*feta.^2);                            % SRC滤波器的调频率，式6.22        
Km = Kr./(1-Kr./Ksrc);                                                      % 改变后的距离向调频率，式6.21
%clear Ksrc
RCMbulk = (1./D_feta_Vref-1/D_fref_Vref)*Rref;                              % 一致RCM，，式7.22
alpha = D_fref_Vref./D_feta_Vref-1;                                         % 频偏参数，式7.28

%% 对齐到方位向多普勒中心
S0 = data.*exp(-1i*2*pi*fc*(eta'*ones(1,Nr)));                              % 搬移至多普勒频率中心
%clear data

%% 变换到距离多普勒域，实现变标相乘
Srd = fftshift(fft(fftshift(S0,1),[],1),1);                                 % 原信号做方位向FFT
%clear S0
tt = 2/C*(R0/D_fref_V+RCMbulk)-2*Rref./(C*D_feta_Vref);                     % P205 (7.26) (7.27)
Ssc = exp(1i*pi*Km.*alpha.*tt.^2);                                          % 变标方程 P207 (7.30)
%clear tt
S1 = Srd.*(Ssc'*ones(1,Nr));                                                % 变标相乘，，式7.31
%clear Srd Ssc

%% 变换到二维频域，实现RC、SRC和一致RCMC
S2 = fftshift(fft(fftshift(S1,2),[],2),2);                                  % 信号变换到二维频域
%clear S1
WindowR = ones(Na,1)*kaiser(Nr,2.5)';                                       % 距离窗
WindowA = kaiser(Na,2.5)*ones(1,Nr);                                        % 方位窗
S2 = S2.*WindowR.*WindowA;                                                  % 加窗，式7.32
%clear WindowR WindowA
Hm = exp(1i*pi./((Km'*ones(1,Nr)).*(1+alpha'*ones(1,Nr))).*(ones(Na,1)*ftau).^2); % 合并的距离压缩和一致RCMC滤波器，，式7.32,1+alpha
%clear alpha
Hrcm = exp(1i*4*pi/C*(RCMbulk'*ones(1,Nr)).*(ones(Na,1)*ftau));             % SRC滤波器,，式7.32,，式7.22
%clear RCMbulk
S3 = S2.*Hm.*Hrcm;                                                          % 相位相乘，，式7.33
%clear S2 Hm Hrcm

%% 变换到距离多普勒域，实现方位压缩和附加相位校正
S4 = ifftshift(ifft(ifftshift(S3,2),[],2),2);                               % 距离向IFFT，式7.34
%clear S3

figure,imagesc(abs(S4));axis image;set(gcf,'Color','w');
title('RD域：RC、SRC和一致RCMC后的信号幅度');
xlabel('距离向（采样点）');ylabel('方位向（采样点）');

Hac = exp(-1i*4*pi*R0*f0*D_feta_V/C);                                       % 方位匹配滤波器，式7.34第一个指数项
Hpc = exp(1i*(4*pi*Km/C^2).*(1-D_feta_Vref/D_fref_Vref).*(R0./D_feta_V-Rref./D_feta_V).^2); % 相位校正滤波器，式7.34第二个指数项
S5 = S4.*(Hac'*ones(1,Nr)).*(Hpc'*ones(1,Nr));                              % 相位相乘，实现滤波，或者取H复共轭相消
%clear S4 Hac Hpc Km

%% 变换到图像域
Stt = ifftshift(ifft(ifftshift(S5,1),[],1),1);                              % 信号的方位向IFFT
%clear S5

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

imwrite(Img,'SAR_CSA.bmp','bmp');
