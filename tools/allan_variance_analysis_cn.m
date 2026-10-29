%% IMU Allan 方差分析
% 读取 output/成果_REAL_松组合真实算法/imu_decoded.txt。
% 输出：
%   output/成果_MATLAB分析/Allan方差分析/allan_gyro.png
%   output/成果_MATLAB分析/Allan方差分析/allan_acc.png
%   output/成果_MATLAB分析/Allan方差分析/allan_summary.txt

clear; clc; close all;

scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    analysisDir = pwd;
else
    analysisDir = fileparts(scriptPath);
end
outputDir = fileparts(analysisDir);
imuFile = fullfile(outputDir, '成果_REAL_松组合真实算法', 'imu_decoded.txt');
outDir = fullfile(analysisDir, 'Allan方差分析');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

checkFile(imuFile);
imu = readImuDecoded(imuFile);
dt = median(diff(imu.tow));
fs = 1 / dt;

gyro = [imu.dthetaX, imu.dthetaY, imu.dthetaZ] ./ dt;
acc = [imu.dvelX, imu.dvelY, imu.dvelZ] ./ dt;

% 为避免动态运动段污染 Allan 曲线，优先使用前 120 s 数据。
% 如果数据长度不足，则使用全部数据。
t0 = imu.tow(1);
staticMask = imu.tow <= t0 + 120;
if sum(staticMask) < 2000
    staticMask = true(size(imu.tow));
end

gyroSeg = gyro(staticMask, :);
accSeg = acc(staticMask, :);
segSeconds = sum(staticMask) * dt;

[tauG, adevG] = allanDeviation(gyroSeg, dt);
[tauA, adevA] = allanDeviation(accSeg, dt);

gyroArw = estimateWhiteNoise(tauG, adevG);
accVrw = estimateWhiteNoise(tauA, adevA);
gyroBiasInstability = min(adevG, [], 1);
accBiasInstability = min(adevA, [], 1);

%% 图 1：陀螺 Allan deviation
f = figure('Name', '陀螺Allan方差', 'Color', 'w');
loglog(tauG, adevG(:, 1), 'r-', tauG, adevG(:, 2), 'g-', tauG, adevG(:, 3), 'b-', 'LineWidth', 1.2);
grid on;
xlabel('平均时间 \tau (s)');
ylabel('Allan deviation (rad/s)');
title('陀螺 Allan deviation');
legend('gyro X', 'gyro Y', 'gyro Z', 'Location', 'best');
exportgraphics(f, fullfile(outDir, 'allan_gyro.png'), 'Resolution', 200);

%% 图 2：加速度计 Allan deviation
f = figure('Name', '加速度计Allan方差', 'Color', 'w');
loglog(tauA, adevA(:, 1), 'r-', tauA, adevA(:, 2), 'g-', tauA, adevA(:, 3), 'b-', 'LineWidth', 1.2);
grid on;
xlabel('平均时间 \tau (s)');
ylabel('Allan deviation (m/s^2)');
title('加速度计 Allan deviation');
legend('acc X', 'acc Y', 'acc Z', 'Location', 'best');
exportgraphics(f, fullfile(outDir, 'allan_acc.png'), 'Resolution', 200);

%% 文本报告
summaryFile = fullfile(outDir, 'allan_summary.txt');
fid = fopen(summaryFile, 'w');
fprintf(fid, 'IMU Allan 方差分析结果\n');
fprintf(fid, '=====================\n\n');
fprintf(fid, '输入文件：%s\n', imuFile);
fprintf(fid, '采样间隔：%.6f s\n', dt);
fprintf(fid, '采样频率：%.3f Hz\n', fs);
fprintf(fid, '用于分析的数据长度：%.3f s\n', segSeconds);
fprintf(fid, '说明：当前使用动态实验开头数据进行 Allan 诊断，不等价于完整 2 小时静态 Allan 标定。\n');
fprintf(fid, '      结果可用于报告中的短时噪声观察；若要严格评估 IMU 随机误差，应补充 2 小时静态数据。\n\n');

fprintf(fid, '陀螺白噪声近似 ARW，单位 rad/s/sqrt(Hz)：\n');
fprintf(fid, '  X %.8e\n  Y %.8e\n  Z %.8e\n\n', gyroArw(1), gyroArw(2), gyroArw(3));
fprintf(fid, '加速度计白噪声近似 VRW，单位 m/s^2/sqrt(Hz)：\n');
fprintf(fid, '  X %.8e\n  Y %.8e\n  Z %.8e\n\n', accVrw(1), accVrw(2), accVrw(3));

fprintf(fid, '陀螺 Allan deviation 最小值，单位 rad/s：\n');
fprintf(fid, '  X %.8e\n  Y %.8e\n  Z %.8e\n\n', gyroBiasInstability(1), gyroBiasInstability(2), gyroBiasInstability(3));
fprintf(fid, '加速度计 Allan deviation 最小值，单位 m/s^2：\n');
fprintf(fid, '  X %.8e\n  Y %.8e\n  Z %.8e\n', accBiasInstability(1), accBiasInstability(2), accBiasInstability(3));
fclose(fid);

fprintf('Allan 方差分析完成。\n');
fprintf('输出目录：%s\n', outDir);

function checkFile(path)
    if ~exist(path, 'file')
        error('找不到文件：%s', path);
    end
end

function imu = readImuDecoded(path)
    fid = fopen(path, 'r');
    cleaner = onCleanup(@() fclose(fid));
    c = textscan(fid, '%f %f %f %f %f %f %f', ...
        'CommentStyle', '#', 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
    imu.tow = c{1};
    imu.dthetaX = c{2};
    imu.dthetaY = c{3};
    imu.dthetaZ = c{4};
    imu.dvelX = c{5};
    imu.dvelY = c{6};
    imu.dvelZ = c{7};
end

function [tau, adev] = allanDeviation(x, dt)
    n = size(x, 1);
    maxM = floor(n / 10);
    mList = unique(floor(logspace(0, log10(maxM), 60)));
    mList = mList(mList >= 1);
    tau = mList(:) * dt;
    adev = nan(numel(mList), size(x, 2));
    for k = 1:numel(mList)
        m = mList(k);
        nb = floor(n / m);
        if nb < 3
            continue;
        end
        y = x(1:nb*m, :);
        y = reshape(y, m, nb, size(x, 2));
        ybar = squeeze(mean(y, 1));
        if nb == 1
            ybar = ybar(:)';
        end
        dy = diff(ybar, 1, 1);
        adev(k, :) = sqrt(0.5 * mean(dy.^2, 1));
    end
end

function noise = estimateWhiteNoise(tau, adev)
    noise = nan(1, size(adev, 2));
    idx = tau >= 0.02 & tau <= 1.0;
    for j = 1:size(adev, 2)
        y = adev(idx, j) .* sqrt(tau(idx));
        y = y(isfinite(y));
        if isempty(y)
            noise(j) = NaN;
        else
            noise(j) = median(y);
        end
    end
end
