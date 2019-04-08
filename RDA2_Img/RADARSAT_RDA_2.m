% RADARSAT-1的数据处理之RDA成像SRC2
% 精确实现二次距离压缩的RDA

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
deta = 1/Fa;                                                                % 方位向采样时间间隔
% 其他参数
t0 = 6.5956e-3;                                                             % 获取数据时间
R0 = t0*C/2;                                                                % 最短斜距

%% 获取原始数据
oldFolder = cd('.\scene01');
load CDdata1.mat;
cd(oldFolder);
a=1;
data=double(data);
[Na,Nr] = size(data);

%% 确定参考距离和参考时间
Rref = (t0+Nr/2*dtau)*C/2;                                                  % 参考距离选在成像区域中心

%% 数据末尾补零
Za = 800;                                                                   % 方位向补零数
Zr = ceil(Tr/dtau);                                                         % 距离向补零数，每个点被波束照射Tr
data = cat(2,data,zeros(Na,Zr));                                            % 距离向补零(右边)
data = cat(1,zeros(Za,Nr+Zr),data);                                         % 方位向补零(上方)
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

%% 对齐到方位向多普勒中心
S0 = data.*exp(-1i*2*pi*fc*(eta'*ones(1,Nr)));                              % 搬移至多普勒频率中心(式5.19)
clear data

%% RC
% 方式2：复制脉冲补零后进行DFT,对结果取复共轭（无时间反褶）
Tref = t0+Nr/2*dtau;                                                        % 距离向参考时间
ht_rc = (abs((ones(Na,1)*tau)-Tref)<Tr/2).*exp(1i*pi*Kr*((ones(Na,1)*tau)-Tref).^2); % 方式2的匹配滤波器
Hf_rc = conj(fftshift(fft(fftshift(ht_rc,2),[],2),2));                      % 滤波器的距离向FFT变换
clear ht_rc
S0_rc = fftshift(fft(fftshift(S0,2),[],2),2);                               % 信号的距离向FFT变换
clear S0
S1 = ifftshift(ifft(ifftshift(S0_rc.*Hf_rc,2),[],2),2);                     % 匹配滤波后距离向IFFT变换
clear Hf_rc S0_rc

figure,imagesc(abs(S1));axis image;set(gcf,'Color','w');
title('时域：距离压缩后的信号幅度');
xlabel('距离向（采样点）');ylabel('方位向（采样点）');

%% SRC二次距离压缩
Sff = fftshift(fft2(fftshift(S1)));                                         % 信号的二维FFT变换
WindowR = ones(Na,1)*kaiser(Nr,2.5)';                                       % 距离窗
WindowA = kaiser(Na,2.5)*ones(1,Nr);                                        % 方位窗
Sff = Sff.*WindowR.*WindowA;                                                % 加窗
clear WindowR WindowA
D = sqrt(1-lambda^2*feta.^2/(4*V^2));                                       % 距离徙动因子，式6.24
Ksrc = (2*V^2*f0^3*(D'*ones(1,Nr)).^3)./(C*R0*(feta'*ones(1,Nr)).^2);       % SRC滤波器的调频率，式6.22
Hsrc = exp(-1i*pi*((ones(Na,1)*ftau).^2)./Ksrc);                            % SRC滤波器，式6.27
clear Ksrc
S1 = Sff.*Hsrc;                                                             % 滤波
clear Sff Hsrc
Srd = ifftshift(ifft(ifftshift(S1,2),[],2),2);                              % SRC校正后距离向IFFT变换，变换到距离多普勒域
clear S1

figure,imagesc(abs(Srd));axis image;set(gcf,'Color','w');
title('RD域：SRC后的信号幅度');
xlabel('距离向（采样点）');ylabel('方位向（采样点）');

%% RCMC
deltaR = dtau*C/2;                                                          % 距离向采样长度间隔
deltaRCM = Rref*(1./D-1);                                                   % 距离徙动量，式6.25
deltaRCM = deltaRCM/deltaR;                                                 % 距离徙动量的量化
IntNum = ceil(deltaRCM);                                                    % 对距离徙动量向上取整;
DecNum = IntNum-deltaRCM;                                                   % 距离徙动量的小数部分量化为（1/16）的倍数
Srcmc = zeros(Na,Nr);
h = waitbar(0,'Sinc插值');
for ii = 1:Na
    SincVal = sinc(DecNum(ii)-4:1:DecNum(ii)+3)';                           % 生成插值核
    for jj = 1:Nr
        kk = jj+IntNum(ii);                                                 % 原矩阵中距离徙动的位置
        if(5<=kk && kk<=Nr-3)                                               % 防止下标溢出
            Srcmc(ii,jj) = Srd(ii,kk-4:kk+3)*SincVal;                       % 插值
        end
    end
    waitbar(ii/Na);
end
close(h);
clear deltaRCM IntNum DecNum Srd

figure,imagesc(abs(Srcmc));axis image;set(gcf,'Color','w');
title('RD域：RCMC后的信号幅度');
xlabel('距离向（采样点）');ylabel('方位向（采样点）');

%% AC
Haz = exp(1i*4*pi*R0*(D'*ones(1,Nr))*f0/C);                                 % 方位向匹配滤波器，式6.26
clear D
S2 = Srcmc.*Haz;                                                            % 匹配滤波
clear Haz Srcmc
Stt = ifftshift(ifft(ifftshift(S2,1),[],1),1);                              % 方位向IFFT变换
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

imwrite(Img,'SAR_RDA_2.bmp','bmp');
