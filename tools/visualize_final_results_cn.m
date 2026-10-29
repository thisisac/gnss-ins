%% 最终成果可视化与统计分析
% 使用方法：
% 1. 在 MATLAB 中打开本文件。
% 2. 直接点击运行。
% 3. 图片会保存到 output/成果_MATLAB分析/最终可视化图片。
% 4. 文字统计会保存到 output/成果_MATLAB分析/最终可视化统计报告.txt。

clear; clc; close all;

%% 路径设置
scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    analysisDir = pwd;
else
    analysisDir = fileparts(scriptPath);
end

outputDir = fileparts(analysisDir);
rootDir = fileparts(outputDir);

rtkFile = fullfile(outputDir, '成果_REAL_松组合真实算法', 'gnss_rtk.pos');
insDir = fullfile(outputDir, '成果_REAL_松组合真实算法');

forwardFile = fullfile(insDir, 'gnss_ins_forward.nav');
forwardSmoothFile = fullfile(insDir, 'gnss_ins_forward_smooth.nav');
backwardFile = fullfile(insDir, 'gnss_ins_backward.nav');
smoothFile = fullfile(insDir, 'gnss_ins_smooth.nav');

refRtkFile = fullfile(rootDir, 'reference', 'rtk', 'gnss_20260604_105437_599404.pos');
refInsFile = fullfile(rootDir, 'reference', 'GNSSINS', 'gnss_20260604_105437_599404.pos');

figDir = fullfile(analysisDir, '最终可视化图片_真实算法重跑');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

reportFile = fullfile(analysisDir, '最终可视化统计报告.txt');

checkFile(rtkFile);
checkFile(forwardFile);
checkFile(forwardSmoothFile);
checkFile(backwardFile);
checkFile(smoothFile);
checkFile(refRtkFile);
checkFile(refInsFile);

%% 读取数据
fprintf('正在读取 RTK 结果...\n');
rtk = readOurRtk(rtkFile);

fprintf('正在读取 GNSS/INS 前向结果...\n');
forward = readOurNav(forwardFile);

fprintf('正在读取 GNSS/INS 前向平滑结果...\n');
forwardSmooth = readOurNav(forwardSmoothFile);

fprintf('正在读取 GNSS/INS 后向结果...\n');
backward = readOurNav(backwardFile);

fprintf('正在读取 GNSS/INS 前后向平滑结果...\n');
smooth = readOurNav(smoothFile);

fprintf('正在读取参考 RTK 与参考 GNSS/INS...\n');
refRtk = readRefPos(refRtkFile, false);
refIns = readRefPos(refInsFile, true);

fprintf('\n数据读取完成：\n');
fprintf('  RTK                 %d 历元\n', numel(rtk.tow));
fprintf('  前向 INS            %d 历元\n', numel(forward.tow));
fprintf('  前向平滑 INS        %d 历元\n', numel(forwardSmooth.tow));
fprintf('  后向 INS            %d 历元\n', numel(backward.tow));
fprintf('  前后向平滑 INS      %d 历元\n', numel(smooth.tow));
fprintf('  参考 RTK            %d 历元\n', numel(refRtk.tow));
fprintf('  参考 GNSS/INS       %d 历元\n', numel(refIns.tow));

%% 误差统计
[rtkM, refRtkM] = matchByTow(rtk, refRtk, 1e-3);
rtkErr = calcErrorEnu(rtkM, refRtkM);

[forwardM, refForwardM] = matchByNearestSecond(forward, refIns, 0.01);
[forwardSmoothM, refForwardSmoothM] = matchByNearestSecond(forwardSmooth, refIns, 0.01);
[backwardM, refBackwardM] = matchByNearestSecond(backward, refIns, 0.01);
[smoothM, refSmoothM] = matchByNearestSecond(smooth, refIns, 0.01);

errForward = calcErrorEnu(forwardM, refForwardM);
errForwardSmooth = calcErrorEnu(forwardSmoothM, refForwardSmoothM);
errBackward = calcErrorEnu(backwardM, refBackwardM);
errSmooth = calcErrorEnu(smoothM, refSmoothM);

statsRtkAll = errorStats(rtkErr);
statsRtkFix = errorStats(subsetError(rtkErr, rtkM.q == 1));
statsRtkFloat = errorStats(subsetError(rtkErr, rtkM.q == 2));

statsForward = errorStats(errForward);
statsForwardSmooth = errorStats(errForwardSmooth);
statsBackward = errorStats(errBackward);
statsSmooth = errorStats(errSmooth);

%% 输出文字报告
fid = fopen(reportFile, 'w');
fprintf(fid, 'GNSS/INS 最终成果可视化统计报告\n');
fprintf(fid, '================================\n\n');
fprintf(fid, '一、读取文件\n');
fprintf(fid, 'RTK 结果：%s\n', rtkFile);
fprintf(fid, '前向结果：%s\n', forwardFile);
fprintf(fid, '前向平滑：%s\n', forwardSmoothFile);
fprintf(fid, '后向结果：%s\n', backwardFile);
fprintf(fid, '前后向平滑：%s\n', smoothFile);
fprintf(fid, 'RTK 参考：%s\n', refRtkFile);
fprintf(fid, '组合导航参考：%s\n\n', refInsFile);
fprintf(fid, '说明：参考软件 pos 文件中的经纬度字段只显示到“度 分”整数，不能用于厘米级水平误差统计。\n');
fprintf(fid, '      因此本报告只用参考结果统计高程误差和速度误差。\n');
fprintf(fid, '      轨迹图使用参考文件中的 Local-VE、Local-VN 速度积分得到相对轨迹，并与本程序轨迹在同一有效历元对齐。\n\n');

fprintf(fid, '二、数据数量\n');
fprintf(fid, 'RTK：%d 历元\n', numel(rtk.tow));
fprintf(fid, '前向 INS：%d 历元\n', numel(forward.tow));
fprintf(fid, '前向平滑 INS：%d 历元\n', numel(forwardSmooth.tow));
fprintf(fid, '后向 INS：%d 历元\n', numel(backward.tow));
fprintf(fid, '前后向平滑 INS：%d 历元\n\n', numel(smooth.tow));

fprintf(fid, '三、RTK 质量统计\n');
fprintf(fid, 'q=1 固定解：%d\n', sum(rtk.q == 1));
fprintf(fid, 'q=2 浮点解：%d\n', sum(rtk.q == 2));
fprintf(fid, 'q=5 fallback：%d\n\n', sum(rtk.q == 5));

fprintf(fid, '四、RTK 对比参考结果\n');
writeErrorStats(fid, 'RTK 全部匹配历元', statsRtkAll);
writeErrorStats(fid, 'RTK 固定解 q=1', statsRtkFix);
writeErrorStats(fid, 'RTK 浮点解 q=2', statsRtkFloat);

fprintf(fid, '\n五、GNSS/INS 对比参考结果\n');
writeErrorStats(fid, '前向松组合', statsForward);
writeErrorStats(fid, '前向平滑', statsForwardSmooth);
writeErrorStats(fid, '后向解算', statsBackward);
writeErrorStats(fid, '前后向平滑', statsSmooth);

fprintf(fid, '\n六、生成图片\n');
fprintf(fid, '01_RTK质量与高程误差.png\n');
fprintf(fid, '02_RTK轨迹与固定解分布.png\n');
fprintf(fid, '03_组合导航轨迹对比.png\n');
fprintf(fid, '04_组合导航高程误差对比.png\n');
fprintf(fid, '05_组合导航速度误差对比.png\n');
fprintf(fid, '06_前后向平滑姿态速度.png\n');
fprintf(fid, '07_三维局部轨迹.png\n');
fprintf(fid, '08_结果与参考直接曲线对比.png\n');
fclose(fid);

fprintf('\n统计报告已保存：%s\n', reportFile);

%% 图 1：RTK 质量与高程误差
f = figure('Name', 'RTK质量与高程误差', 'Color', 'w');
t = rtkM.tow - rtkM.tow(1);
tiledlayout(3, 1, 'TileSpacing', 'compact');

nexttile;
stairs(t, rtkM.q, 'LineWidth', 1.1);
grid on; ylim([0.5 5.5]);
ylabel('质量标记 q');
title('RTK 解算质量');

nexttile;
plot(t, rtkM.ns, '.-', 'MarkerSize', 6);
grid on;
ylabel('卫星数');

nexttile; hold on;
idxFix = rtkM.q == 1;
idxFloat = rtkM.q == 2;
idxFallback = rtkM.q == 5;
plot(t(idxFallback), rtkErr.dU(idxFallback), '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', 6);
plot(t(idxFloat), rtkErr.dU(idxFloat), '.', 'Color', [0.85 0.45 0.15], 'MarkerSize', 7);
plot(t(idxFix), rtkErr.dU(idxFix), 'b.', 'MarkerSize', 8);
yline(0, 'k-');
grid on;
xlabel('相对时间 (s)');
ylabel('高程误差 (m)');
legend('fallback', 'float', 'fixed', 'Location', 'best');
saveFigure(f, figDir, '01_RTK质量与高程误差.png');

%% 图 2：RTK 轨迹与参考轨迹对比
anchorRtk = find(rtkM.q == 1, 1, 'first');
if isempty(anchorRtk)
    anchorRtk = 1;
end
[eRtk, nRtk, ~] = llhToLocal(rtkM.lat, rtkM.lon, rtkM.h, rtkM.lat(anchorRtk), rtkM.lon(anchorRtk), rtkM.h(anchorRtk));
refRtkTrack = integrateVelocityTrack(refRtkM, anchorRtk);

f = figure('Name', 'RTK轨迹与参考轨迹对比', 'Color', 'w');
plot(refRtkTrack.e, refRtkTrack.n, 'k-', 'LineWidth', 1.2); hold on;
plot(eRtk(rtkM.q ~= 1), nRtk(rtkM.q ~= 1), '.', 'Color', [0.70 0.70 0.70], 'MarkerSize', 5);
plot(eRtk(rtkM.q == 1), nRtk(rtkM.q == 1), 'r.', 'MarkerSize', 9);
axis equal; grid on;
xlabel('东向 (m)');
ylabel('北向 (m)');
title('RTK 轨迹与参考轨迹对比');
legend('参考 RTK 速度积分轨迹', '自研非固定解', '自研固定解', 'Location', 'best');
saveFigure(f, figDir, '02_RTK轨迹与固定解分布.png');

%% 图 3：组合导航轨迹与参考轨迹对比
anchorIns = 1;
[eF, nF, ~] = llhToLocal(forwardM.lat, forwardM.lon, forwardM.h, smoothM.lat(anchorIns), smoothM.lon(anchorIns), smoothM.h(anchorIns));
[eFS, nFS, ~] = llhToLocal(forwardSmoothM.lat, forwardSmoothM.lon, forwardSmoothM.h, smoothM.lat(anchorIns), smoothM.lon(anchorIns), smoothM.h(anchorIns));
[eB, nB, ~] = llhToLocal(backwardM.lat, backwardM.lon, backwardM.h, smoothM.lat(anchorIns), smoothM.lon(anchorIns), smoothM.h(anchorIns));
[eS, nS, ~] = llhToLocal(smoothM.lat, smoothM.lon, smoothM.h, smoothM.lat(anchorIns), smoothM.lon(anchorIns), smoothM.h(anchorIns));
uS = smoothM.h - smoothM.h(anchorIns);
refInsTrack = integrateVelocityTrack(refSmoothM, anchorIns);

f = figure('Name', '组合导航轨迹与参考轨迹对比', 'Color', 'w');
plot(refInsTrack.e, refInsTrack.n, 'k-', 'LineWidth', 1.4); hold on;
plot(eF, nF, 'Color', [0.80 0.80 0.80], 'LineWidth', 0.7);
plot(eFS, nFS, 'Color', [0.25 0.60 0.95], 'LineWidth', 0.9);
plot(eB, nB, 'Color', [0.25 0.70 0.25], 'LineWidth', 0.9);
plot(eS, nS, 'r-', 'LineWidth', 1.2);
markStartEnd(eS, nS);
axis equal; grid on;
xlabel('东向 (m)');
ylabel('北向 (m)');
title('真实算法 GNSS/INS 松组合轨迹与参考轨迹对比');
legend('参考 GNSS/INS 速度积分轨迹', '前向', '前向平滑', '后向', '前后向平滑', '起点', '终点', 'Location', 'best');
saveFigure(f, figDir, '03_组合导航轨迹对比.png');

%% 图 4：组合导航高程误差对比
f = figure('Name', '组合导航高程误差对比', 'Color', 'w');
plot(forwardM.tow - refIns.tow(1), errForward.dU, '-', 'Color', [0.65 0.65 0.65]); hold on;
plot(forwardSmoothM.tow - refIns.tow(1), errForwardSmooth.dU, '-', 'Color', [0.25 0.60 0.95]);
plot(backwardM.tow - refIns.tow(1), errBackward.dU, '-', 'Color', [0.25 0.70 0.25]);
plot(smoothM.tow - refIns.tow(1), errSmooth.dU, 'r-', 'LineWidth', 1.2);
yline(0, 'k-');
grid on;
xlabel('相对参考起点时间 (s)');
ylabel('高程误差 (m)');
title('组合导航高程误差');
legend('前向', '前向平滑', '后向', '前后向平滑', 'Location', 'best');
saveFigure(f, figDir, '04_组合导航高程误差对比.png');

%% 图 5：组合导航速度误差对比
f = figure('Name', '组合导航速度误差对比', 'Color', 'w');
plot(forwardM.tow - refIns.tow(1), errForward.v3, '-', 'Color', [0.65 0.65 0.65]); hold on;
plot(forwardSmoothM.tow - refIns.tow(1), errForwardSmooth.v3, '-', 'Color', [0.25 0.60 0.95]);
plot(backwardM.tow - refIns.tow(1), errBackward.v3, '-', 'Color', [0.25 0.70 0.25]);
plot(smoothM.tow - refIns.tow(1), errSmooth.v3, 'r-', 'LineWidth', 1.2);
grid on;
xlabel('相对参考起点时间 (s)');
ylabel('三维速度误差 (m/s)');
title('组合导航三维速度误差');
legend('前向', '前向平滑', '后向', '前后向平滑', 'Location', 'best');
saveFigure(f, figDir, '05_组合导航速度误差对比.png');

%% 图 6：前后向平滑姿态速度与参考对比
f = figure('Name', '前后向平滑姿态速度与参考对比', 'Color', 'w');
tSmooth = smoothM.tow - smoothM.tow(1);
tRef = refSmoothM.tow - refSmoothM.tow(1);
tiledlayout(2, 1, 'TileSpacing', 'compact');

nexttile;
plot(tSmooth, smoothM.roll, 'b-', tSmooth, smoothM.pitch, 'r-', tSmooth, wrapTo180Local(smoothM.yaw), 'Color', [0.95 0.55 0.05]); hold on;
plot(tRef, refSmoothM.roll, 'b--', tRef, refSmoothM.pitch, 'r--', tRef, wrapTo180Local(refSmoothM.yaw), '--', 'Color', [0.45 0.45 0.45]);
grid on;
ylabel('姿态角 (deg)');
title('前后向平滑姿态与参考对比');
legend('roll', 'pitch', 'yaw', '参考 roll', '参考 pitch', '参考 yaw', 'Location', 'best');

nexttile;
speedSmooth = sqrt(smoothM.vn.^2 + smoothM.ve.^2 + smoothM.vd.^2);
speedRef = sqrt(refSmoothM.vn.^2 + refSmoothM.ve.^2 + refSmoothM.vd.^2);
plot(tSmooth, speedSmooth, 'r-', 'LineWidth', 1.2); hold on;
plot(tRef, speedRef, 'k--', 'LineWidth', 1.1);
grid on;
xlabel('相对时间 (s)');
ylabel('速度 (m/s)');
title('前后向平滑速度模长与参考对比');
legend('前后向平滑', '参考 GNSS/INS', 'Location', 'best');
saveFigure(f, figDir, '06_前后向平滑姿态速度.png');

%% 图 7：三维局部轨迹与参考对比
idx3d = thinIndex(numel(eS), 100000);
f = figure('Name', '三维局部轨迹与参考对比', 'Color', 'w');
plot3(refInsTrack.e, refInsTrack.n, refInsTrack.u, 'k-', 'LineWidth', 1.1); hold on;
plot3(eS(idx3d), nS(idx3d), uS(idx3d), 'r-', 'LineWidth', 1.2);
grid on;
xlabel('东向 (m)');
ylabel('北向 (m)');
zlabel('天向 (m)');
title('真实算法前后向平滑三维轨迹与参考对比');
legend('参考 GNSS/INS 速度积分轨迹', '前后向平滑', 'Location', 'best');
view(35, 25);
saveFigure(f, figDir, '07_三维局部轨迹.png');

%% 图 8：结果与参考直接曲线对比
f = figure('Name', '结果与参考直接曲线对比', 'Color', 'w');
tiledlayout(3, 1, 'TileSpacing', 'compact');

nexttile;
plot(tSmooth, smoothM.h, 'r-', 'LineWidth', 1.1); hold on;
plot(tRef, refSmoothM.h, 'k--', 'LineWidth', 1.0);
grid on;
ylabel('高程 (m)');
title('前后向平滑结果与参考结果直接对比');
legend('前后向平滑', '参考 GNSS/INS', 'Location', 'best');

nexttile;
plot(tSmooth, smoothM.vn, 'b-', tSmooth, smoothM.ve, 'r-', tSmooth, smoothM.vd, 'Color', [0.95 0.55 0.05]); hold on;
plot(tRef, refSmoothM.vn, 'b--', tRef, refSmoothM.ve, 'r--', tRef, refSmoothM.vd, '--', 'Color', [0.45 0.45 0.45]);
grid on;
ylabel('速度分量 (m/s)');
legend('v_N', 'v_E', 'v_D', '参考 v_N', '参考 v_E', '参考 v_D', 'Location', 'best');

nexttile;
plot(tSmooth, smoothM.h - refSmoothM.h, 'r-', 'LineWidth', 1.1);
yline(0, 'k-');
grid on;
xlabel('相对时间 (s)');
ylabel('高程差 (m)');
saveFigure(f, figDir, '08_结果与参考直接曲线对比.png');

fprintf('图片已保存到：%s\n', figDir);
fprintf('全部完成。\n');

%% 本文件使用的局部函数
function checkFile(path)
    if ~exist(path, 'file')
        error('找不到文件：%s', path);
    end
end

function saveFigure(figHandle, figDir, fileName)
    exportgraphics(figHandle, fullfile(figDir, fileName), 'Resolution', 200);
end

function s = readOurRtk(path)
    fid = fopen(path, 'r');
    cleaner = onCleanup(@() fclose(fid));
    fmt = ['%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s'];
    c = textscan(fid, fmt, 'CommentStyle', '#', 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
    s.week = c{1};
    s.tow = c{2};
    s.lat = c{3};
    s.lon = c{4};
    s.h = c{5};
    s.bx = c{6};
    s.by = c{7};
    s.bz = c{8};
    s.vn = c{9};
    s.ve = c{10};
    s.vd = c{11};
    s.q = c{12};
    s.ns = c{13};
    s.ratio = c{14};
    s.stdN = c{15};
    s.stdE = c{16};
    s.stdD = c{17};
    s.mode = string(c{18});
end

function s = readOurNav(path)
    fid = fopen(path, 'r');
    cleaner = onCleanup(@() fclose(fid));
    fmt = ['%f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %*[^\n]'];
    c = textscan(fid, fmt, 'CommentStyle', '#', 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
    s.tow = c{1};
    s.x = c{2};
    s.y = c{3};
    s.z = c{4};
    s.lat = c{5};
    s.lon = c{6};
    s.h = c{7};
    s.vn = c{8};
    s.ve = c{9};
    s.vd = c{10};
    s.roll = c{11};
    s.pitch = c{12};
    s.yaw = c{13};
    s.q = c{14};
    s.mode = string(c{15});
end

function s = readRefPos(path, hasAttitude)
    lines = readlines(path);
    vals = [];
    for i = 1:numel(lines)
        line = strtrim(lines(i));
        if line == ""
            continue;
        end
        firstChar = extractBefore(line + " ", 2);
        if ~isstrprop(firstChar, 'digit')
            continue;
        end
        p = split(line);
        if numel(p) < 11
            continue;
        end
        tow = str2double(p(3));
        lat = degmin(str2double(p(4)), str2double(p(5)));
        lon = degmin(str2double(p(6)), str2double(p(7)));
        h = str2double(p(8));
        ve = str2double(p(9));
        vn = str2double(p(10));
        vd = -str2double(p(11));
        roll = NaN; pitch = NaN; yaw = NaN;
        if hasAttitude && numel(p) >= 17
            roll = degmin(str2double(p(12)), str2double(p(13)));
            pitch = degmin(str2double(p(14)), str2double(p(15)));
            yaw = degmin(str2double(p(16)), str2double(p(17)));
        end
        if all(isfinite([tow, lat, lon, h]))
            vals(end + 1, :) = [tow, lat, lon, h, vn, ve, vd, roll, pitch, yaw]; %#ok<AGROW>
        end
    end
    s.tow = vals(:, 1);
    s.lat = vals(:, 2);
    s.lon = vals(:, 3);
    s.h = vals(:, 4);
    s.vn = vals(:, 5);
    s.ve = vals(:, 6);
    s.vd = vals(:, 7);
    s.roll = vals(:, 8);
    s.pitch = vals(:, 9);
    s.yaw = vals(:, 10);
end

function d = degmin(degValue, minuteValue)
    if degValue < 0
        d = degValue - minuteValue / 60;
    else
        d = degValue + minuteValue / 60;
    end
end

function [a2, b2] = matchByTow(a, b, tol)
    ia = [];
    ib = [];
    j = 1;
    for i = 1:numel(a.tow)
        while j < numel(b.tow) && b.tow(j) < a.tow(i) - tol
            j = j + 1;
        end
        if j <= numel(b.tow) && abs(b.tow(j) - a.tow(i)) <= tol
            ia(end + 1, 1) = i; %#ok<AGROW>
            ib(end + 1, 1) = j; %#ok<AGROW>
        end
    end
    a2 = subsetStruct(a, ia);
    b2 = subsetStruct(b, ib);
end

function [a2, b2] = matchByNearestSecond(a, b, tol)
    ia = zeros(numel(b.tow), 1);
    ib = zeros(numel(b.tow), 1);
    n = 0;
    for j = 1:numel(b.tow)
        [d, i] = min(abs(a.tow - b.tow(j)));
        if d <= tol
            n = n + 1;
            ia(n) = i;
            ib(n) = j;
        end
    end
    ia = ia(1:n);
    ib = ib(1:n);
    a2 = subsetStruct(a, ia);
    b2 = subsetStruct(b, ib);
end

function out = subsetStruct(s, idx)
    names = fieldnames(s);
    for k = 1:numel(names)
        v = s.(names{k});
        if size(v, 1) == numel(s.tow)
            out.(names{k}) = v(idx, :);
        else
            out.(names{k}) = v;
        end
    end
end

function err = calcErrorEnu(a, b)
    [ea, na, ua] = llhToLocal(a.lat, a.lon, a.h, b.lat(1), b.lon(1), b.h(1));
    [eb, nb, ub] = llhToLocal(b.lat, b.lon, b.h, b.lat(1), b.lon(1), b.h(1));
    err.dE = ea - eb;
    err.dN = na - nb;
    err.dU = a.h - b.h;
    err.hor = hypot(err.dE, err.dN);
    err.pos3 = sqrt(err.dE.^2 + err.dN.^2 + (ua - ub).^2);
    err.dVn = a.vn - b.vn;
    err.dVe = a.ve - b.ve;
    err.dVd = a.vd - b.vd;
    err.v3 = sqrt(err.dVn.^2 + err.dVe.^2 + err.dVd.^2);
end

function e2 = subsetError(e, idx)
    names = fieldnames(e);
    for k = 1:numel(names)
        e2.(names{k}) = e.(names{k})(idx);
    end
end

function st = errorStats(e)
    st.n = numel(e.dU);
    st.dE = scalarStats(e.dE);
    st.dN = scalarStats(e.dN);
    st.dU = scalarStats(e.dU);
    st.hor = scalarStats(e.hor);
    st.pos3 = scalarStats(e.pos3);
    st.v3 = scalarStats(e.v3);
end

function st = scalarStats(x)
    x = x(:);
    x = x(isfinite(x));
    st.n = numel(x);
    if isempty(x)
        st.mean = NaN;
        st.rms = NaN;
        st.maxabs = NaN;
    else
        st.mean = mean(x);
        st.rms = sqrt(mean(x.^2));
        st.maxabs = max(abs(x));
    end
end

function writeErrorStats(fid, name, st)
    fprintf(fid, '\n%s：匹配历元 %d\n', name, st.n);
    fprintf(fid, '  高程误差 RMS：%.4f m，最大绝对值：%.4f m\n', st.dU.rms, st.dU.maxabs);
    fprintf(fid, '  三维速度 RMS：%.4f m/s，最大值：%.4f m/s\n', st.v3.rms, st.v3.maxabs);
end

function [e, n, u] = llhToLocal(latDeg, lonDeg, h, lat0Deg, lon0Deg, h0)
    r0 = llh2ecefLocal(lat0Deg, lon0Deg, h0);
    r = llh2ecefLocal(latDeg, lonDeg, h);
    d = r - r0;
    lat0 = deg2rad(lat0Deg);
    lon0 = deg2rad(lon0Deg);
    c = [-sin(lon0),            cos(lon0),           0; ...
         -sin(lat0)*cos(lon0), -sin(lat0)*sin(lon0), cos(lat0); ...
          cos(lat0)*cos(lon0),  cos(lat0)*sin(lon0), sin(lat0)];
    enu = (c * d')';
    e = enu(:, 1);
    n = enu(:, 2);
    u = enu(:, 3);
end

function r = llh2ecefLocal(latDeg, lonDeg, h)
    a = 6378137.0;
    f = 1 / 298.257223563;
    e2 = f * (2 - f);
    lat = deg2rad(latDeg);
    lon = deg2rad(lonDeg);
    n = a ./ sqrt(1 - e2 .* sin(lat).^2);
    x = (n + h) .* cos(lat) .* cos(lon);
    y = (n + h) .* cos(lat) .* sin(lon);
    z = (n .* (1 - e2) + h) .* sin(lat);
    r = [x(:), y(:), z(:)];
end

function idx = thinIndex(n, maxCount)
    if n <= maxCount
        idx = (1:n)';
    else
        idx = unique(round(linspace(1, n, maxCount)))';
    end
end

function track = integrateVelocityTrack(s, anchorIdx)
    if nargin < 2 || isempty(anchorIdx)
        anchorIdx = 1;
    end
    tow = s.tow(:);
    ve = s.ve(:);
    vn = s.vn(:);
    vd = s.vd(:);
    n = numel(tow);
    e = zeros(n, 1);
    north = zeros(n, 1);
    u = zeros(n, 1);
    for i = 2:n
        dt = tow(i) - tow(i - 1);
        if ~isfinite(dt) || dt <= 0 || dt > 10
            dt = 0;
        end
        e(i) = e(i - 1) + 0.5 * (ve(i - 1) + ve(i)) * dt;
        north(i) = north(i - 1) + 0.5 * (vn(i - 1) + vn(i)) * dt;
        u(i) = u(i - 1) - 0.5 * (vd(i - 1) + vd(i)) * dt;
    end
    anchorIdx = max(1, min(anchorIdx, n));
    track.e = e - e(anchorIdx);
    track.n = north - north(anchorIdx);
    track.u = s.h(:) - s.h(anchorIdx);
end

function y = wrapTo180Local(x)
    y = mod(x + 180, 360) - 180;
end

function markStartEnd(e, n)
    plot(e(1), n(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 7);
    plot(e(end), n(end), 'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 7);
end
