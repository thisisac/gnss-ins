%% 纯惯导结果分析
% 读取 output/成果_PureINS验证/ins_only.nav，与 reference/GNSSINS 参考结果对比。
% 输出：
%   output/成果_MATLAB分析/纯惯导分析/pure_ins_error.png
%   output/成果_MATLAB分析/纯惯导分析/pure_ins_trajectory.png
%   output/成果_MATLAB分析/纯惯导分析/pure_ins_summary.txt

clear; clc; close all;

scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    analysisDir = pwd;
else
    analysisDir = fileparts(scriptPath);
end
outputDir = fileparts(analysisDir);
rootDir = fileparts(outputDir);

insFile = fullfile(outputDir, '成果_PureINS验证', 'ins_only.nav');
refFile = fullfile(rootDir, 'reference', 'GNSSINS', 'gnss_20260604_105437_599404.pos');
outDir = fullfile(analysisDir, '纯惯导分析');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

checkFile(insFile);
checkFile(refFile);

ins = readOurNav(insFile);
ref = readRefPos(refFile);
[insM, refM] = matchByNearestSecond(ins, ref, 0.01);

[eIns, nIns, uIns] = llhToLocal(insM.lat, insM.lon, insM.h, insM.lat(1), insM.lon(1), insM.h(1));
refTrack = integrateVelocityTrack(refM, 1);

dE = eIns - refTrack.e;
dN = nIns - refTrack.n;
dU = insM.h - refM.h;
hor = hypot(dE, dN);
pos3 = sqrt(dE.^2 + dN.^2 + dU.^2);
v3 = sqrt((insM.vn - refM.vn).^2 + (insM.ve - refM.ve).^2 + (insM.vd - refM.vd).^2);
t = insM.tow - insM.tow(1);

stats = makeStats(dU, hor, pos3, v3);

%% 图 1：误差发散
f = figure('Name', '纯惯导误差发散', 'Color', 'w');
tiledlayout(3, 1, 'TileSpacing', 'compact');

nexttile;
plot(t, dU, 'r-', 'LineWidth', 1.0);
grid on;
ylabel('高程误差 (m)');
title('纯惯导相对参考结果的误差发散');

nexttile;
plot(t, hor, 'b-', 'LineWidth', 1.0);
grid on;
ylabel('水平误差 (m)');

nexttile;
plot(t, v3, 'k-', 'LineWidth', 1.0);
grid on;
xlabel('相对时间 (s)');
ylabel('三维速度误差 (m/s)');
exportgraphics(f, fullfile(outDir, 'pure_ins_error.png'), 'Resolution', 200);

%% 图 2：轨迹对比
f = figure('Name', '纯惯导轨迹对比', 'Color', 'w');
plot(refTrack.e, refTrack.n, 'k-', 'LineWidth', 1.2); hold on;
plot(eIns, nIns, 'r-', 'LineWidth', 1.0);
axis equal; grid on;
xlabel('东向 (m)');
ylabel('北向 (m)');
title('纯惯导轨迹与参考轨迹对比');
legend('参考速度积分轨迹', '纯惯导轨迹', 'Location', 'best');
exportgraphics(f, fullfile(outDir, 'pure_ins_trajectory.png'), 'Resolution', 200);

%% 文本报告
summaryFile = fullfile(outDir, 'pure_ins_summary.txt');
fid = fopen(summaryFile, 'w');
fprintf(fid, '纯惯导机械编排分析结果\n');
fprintf(fid, '=====================\n\n');
fprintf(fid, '输入文件：%s\n', insFile);
fprintf(fid, '参考文件：%s\n', refFile);
fprintf(fid, '匹配历元：%d\n', numel(t));
fprintf(fid, '时长：%.3f s\n\n', t(end) - t(1));
fprintf(fid, '高程误差 RMS：%.4f m，最大绝对值：%.4f m\n', stats.hRms, stats.hMax);
fprintf(fid, '水平误差 RMS：%.4f m，最大值：%.4f m\n', stats.horRms, stats.horMax);
fprintf(fid, '三维位置误差 RMS：%.4f m，最大值：%.4f m\n', stats.pos3Rms, stats.pos3Max);
fprintf(fid, '三维速度误差 RMS：%.4f m/s，最大值：%.4f m/s\n\n', stats.v3Rms, stats.v3Max);
fprintf(fid, '说明：纯惯导没有 GNSS 量测更新，误差会随时间快速发散。\n');
fprintf(fid, '      该结果用于展示机械编排可运行以及 IMU 误差对导航结果的影响。\n');
fclose(fid);

fprintf('纯惯导分析完成。\n');
fprintf('输出目录：%s\n', outDir);

function checkFile(path)
    if ~exist(path, 'file')
        error('找不到文件：%s', path);
    end
end

function s = readOurNav(path)
    fid = fopen(path, 'r');
    cleaner = onCleanup(@() fclose(fid));
    fmt = ['%f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %*[^\n]'];
    c = textscan(fid, fmt, 'CommentStyle', '#', 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
    s.tow = c{1};
    s.lat = c{5};
    s.lon = c{6};
    s.h = c{7};
    s.vn = c{8};
    s.ve = c{9};
    s.vd = c{10};
    s.roll = c{11};
    s.pitch = c{12};
    s.yaw = c{13};
end

function s = readRefPos(path)
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
        if numel(p) < 17
            continue;
        end
        tow = str2double(p(3));
        h = str2double(p(8));
        ve = str2double(p(9));
        vn = str2double(p(10));
        vd = -str2double(p(11));
        if all(isfinite([tow, h, vn, ve, vd]))
            vals(end + 1, :) = [tow, h, vn, ve, vd]; %#ok<AGROW>
        end
    end
    s.tow = vals(:, 1);
    s.h = vals(:, 2);
    s.vn = vals(:, 3);
    s.ve = vals(:, 4);
    s.vd = vals(:, 5);
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

function track = integrateVelocityTrack(s, anchorIdx)
    tow = s.tow(:);
    ve = s.ve(:);
    vn = s.vn(:);
    n = numel(tow);
    e = zeros(n, 1);
    north = zeros(n, 1);
    for i = 2:n
        dt = tow(i) - tow(i - 1);
        if ~isfinite(dt) || dt <= 0 || dt > 10
            dt = 0;
        end
        e(i) = e(i - 1) + 0.5 * (ve(i - 1) + ve(i)) * dt;
        north(i) = north(i - 1) + 0.5 * (vn(i - 1) + vn(i)) * dt;
    end
    anchorIdx = max(1, min(anchorIdx, n));
    track.e = e - e(anchorIdx);
    track.n = north - north(anchorIdx);
end

function st = makeStats(dU, hor, pos3, v3)
    st.hRms = sqrt(mean(dU.^2, 'omitnan'));
    st.hMax = max(abs(dU), [], 'omitnan');
    st.horRms = sqrt(mean(hor.^2, 'omitnan'));
    st.horMax = max(hor, [], 'omitnan');
    st.pos3Rms = sqrt(mean(pos3.^2, 'omitnan'));
    st.pos3Max = max(pos3, [], 'omitnan');
    st.v3Rms = sqrt(mean(v3.^2, 'omitnan'));
    st.v3Max = max(v3, [], 'omitnan');
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
    rn = a ./ sqrt(1 - e2 .* sin(lat).^2);
    x = (rn + h) .* cos(lat) .* cos(lon);
    y = (rn + h) .* cos(lat) .* sin(lon);
    z = (rn .* (1 - e2) + h) .* sin(lat);
    r = [x(:), y(:), z(:)];
end
