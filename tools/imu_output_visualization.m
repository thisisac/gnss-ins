%% IMU output visualization
% Read decoded IMU output and generate report-ready plots.
%
% Default input:
%   output/imu_decoded.txt
%
% Optional input:
%   output/成果_MATLAB分析/IMU活动与对准分析/dynamic_imu_decoded_activity.csv
%
% Outputs:
%   output/成果_MATLAB分析/IMU输出可视化/imu_output_overview.png
%   output/成果_MATLAB分析/IMU输出可视化/imu_output_components.png
%   output/成果_MATLAB分析/IMU输出可视化/imu_output_summary.txt

clear; clc; close all;

scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    analysisDir = pwd;
else
    analysisDir = fileparts(scriptPath);
end
outputDir = fileparts(analysisDir);
rootDir = fileparts(outputDir);

txtFile = fullfile(outputDir, 'imu_decoded.txt');
csvFile = fullfile(outputDir, '成果_MATLAB分析', 'IMU活动与对准分析', 'dynamic_imu_decoded_activity.csv');
alignFile = fullfile(outputDir, 'coarse_alignment_result.txt');
outDir = fullfile(analysisDir, 'IMU输出可视化');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

if exist(txtFile, 'file')
    data = readDecodedImuTxt(txtFile);
    sourceFile = txtFile;
elseif exist(csvFile, 'file')
    data = readDecodedImuCsv(csvFile);
    sourceFile = csvFile;
else
    error('Cannot find IMU decoded file. Expected: %s or %s', txtFile, csvFile);
end

if isempty(data.tow)
    error('IMU data is empty: %s', sourceFile);
end

align = readCoarseAlignmentMeta(alignFile);

t = data.tow - data.tow(1);
dt = data.dt;
if isempty(dt)
    dt = [NaN; diff(data.tow)];
end
dtRate = dt;
validDt = dtRate(isfinite(dtRate) & dtRate > 0);
if isempty(validDt)
    validDt = 1;
end
dtRate(~isfinite(dtRate) | dtRate <= 0) = median(validDt);

% For delta-angle / delta-velocity outputs, derive approximate rates.
gyroRate = data.dtheta ./ dtRate;
accRate = data.dvel ./ dtRate;
gyroNorm = vecnorm(data.dtheta, 2, 2);
accNorm = vecnorm(data.dvel, 2, 2);

%% Overview plot
f1 = figure('Name', 'IMU Output Overview', 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(t, data.dtheta(:,1), 'r-', 'LineWidth', 0.9); hold on;
plot(t, data.dtheta(:,2), 'g-', 'LineWidth', 0.9);
plot(t, data.dtheta(:,3), 'b-', 'LineWidth', 0.9);
if align.hasInfo
    shadePhase(align.alignSeconds, [0.90 0.95 1.00], '静止对准段');
    xline(align.alignSeconds, '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.0);
    text(0.01, 0.93, '蓝色阴影: 静止对准段', 'Units', 'normalized', 'Color', [0.10 0.20 0.55], 'FontSize', 9);
    text(0.01, 0.84, '虚线右侧: 运动段', 'Units', 'normalized', 'Color', [0.20 0.20 0.20], 'FontSize', 9);
end
grid on;
ylabel('d\theta (rad)');
title('IMU delta angle');
legend({'x', 'y', 'z'}, 'Location', 'best');

nexttile;
plot(t, data.dvel(:,1), 'r-', 'LineWidth', 0.9); hold on;
plot(t, data.dvel(:,2), 'g-', 'LineWidth', 0.9);
plot(t, data.dvel(:,3), 'b-', 'LineWidth', 0.9);
if align.hasInfo
    shadePhase(align.alignSeconds, [0.90 0.95 1.00], '');
    xline(align.alignSeconds, '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.0);
end
grid on;
xlabel('Time from start (s)');
ylabel('dV (m/s)');
title('IMU delta velocity');
legend({'x', 'y', 'z'}, 'Location', 'best');

exportgraphics(f1, fullfile(outDir, 'imu_output_overview.png'), 'Resolution', 220);

%% Component plot
f2 = figure('Name', 'IMU Output Components', 'Color', 'w');
tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(t, gyroRate(:,1), 'r-', 'LineWidth', 0.9); hold on;
plot(t, gyroRate(:,2), 'g-', 'LineWidth', 0.9);
plot(t, gyroRate(:,3), 'b-', 'LineWidth', 0.9);
if align.hasInfo
    shadePhase(align.alignSeconds, [0.90 0.95 1.00], '');
    xline(align.alignSeconds, '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.0);
end
grid on;
ylabel('d\theta / dt (rad/s)');
title('Approximate angular rate from delta angle');
legend({'x', 'y', 'z'}, 'Location', 'best');

nexttile;
plot(t, accRate(:,1), 'r-', 'LineWidth', 0.9); hold on;
plot(t, accRate(:,2), 'g-', 'LineWidth', 0.9);
plot(t, accRate(:,3), 'b-', 'LineWidth', 0.9);
if align.hasInfo
    shadePhase(align.alignSeconds, [0.90 0.95 1.00], '');
    xline(align.alignSeconds, '--', 'Color', [0.15 0.15 0.15], 'LineWidth', 1.0);
end
grid on;
ylabel('dV / dt (m/s^2)');
title('Approximate specific force from delta velocity');
legend({'x', 'y', 'z'}, 'Location', 'best');

nexttile;
plot(t, gyroNorm, 'k-', 'LineWidth', 0.9); hold on;
grid on;
ylabel('||d\theta|| (rad)');
title('Delta angle norm');

nexttile;
plot(t, accNorm, 'k-', 'LineWidth', 0.9); hold on;
grid on;
xlabel('Time from start (s)');
ylabel('||dV|| (m/s)');
title('Delta velocity norm');

exportgraphics(f2, fullfile(outDir, 'imu_output_components.png'), 'Resolution', 220);

%% Summary file
summaryFile = fullfile(outDir, 'imu_output_summary.txt');
fid = fopen(summaryFile, 'w');
if fid < 0
    error('Failed to open summary file for writing: %s', summaryFile);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'IMU output visualization summary\n');
fprintf(fid, '===============================\n\n');
fprintf(fid, 'Source file: %s\n', sourceFile);
if align.hasInfo
    fprintf(fid, 'Alignment file: %s\n', alignFile);
    fprintf(fid, 'Static alignment duration: %.3f s\n', align.alignSeconds);
    fprintf(fid, 'Static alignment tow start: %.6f\n', align.imuStartTow);
    fprintf(fid, 'Static alignment tow end: %.6f\n\n', align.imuStartTow + align.alignSeconds);
else
    fprintf(fid, 'Alignment file: not found\n\n');
end
fprintf(fid, 'Samples: %d\n', numel(data.tow));
fprintf(fid, 'Duration: %.3f s\n', t(end));
fprintf(fid, 'Mean dt: %.6f s\n', mean(dt(isfinite(dt))));
fprintf(fid, 'Std dt:  %.6f s\n\n', std(dt(isfinite(dt))));

fprintf(fid, 'Delta angle norm: mean=%.6e rad, std=%.6e rad, max=%.6e rad\n', ...
    mean(gyroNorm), std(gyroNorm), max(gyroNorm));
fprintf(fid, 'Delta velocity norm: mean=%.6e m/s, std=%.6e m/s, max=%.6e m/s\n\n', ...
    mean(accNorm), std(accNorm), max(accNorm));

fprintf(fid, 'Approximate angular rate ranges (rad/s):\n');
fprintf(fid, '  x: [%.6e, %.6e]\n', min(gyroRate(:,1)), max(gyroRate(:,1)));
fprintf(fid, '  y: [%.6e, %.6e]\n', min(gyroRate(:,2)), max(gyroRate(:,2)));
fprintf(fid, '  z: [%.6e, %.6e]\n\n', min(gyroRate(:,3)), max(gyroRate(:,3)));

fprintf(fid, 'Approximate specific force ranges (m/s^2):\n');
fprintf(fid, '  x: [%.6e, %.6e]\n', min(accRate(:,1)), max(accRate(:,1)));
fprintf(fid, '  y: [%.6e, %.6e]\n', min(accRate(:,2)), max(accRate(:,2)));
fprintf(fid, '  z: [%.6e, %.6e]\n', min(accRate(:,3)), max(accRate(:,3)));

fprintf('IMU visualization complete.\n');
fprintf('Source: %s\n', sourceFile);
fprintf('Output directory: %s\n', outDir);

function data = readDecodedImuTxt(path)
    fid = fopen(path, 'r');
    if fid < 0
        error('Cannot open file: %s', path);
    end
    cleanup = onCleanup(@() fclose(fid));
    c = textscan(fid, '%f %f %f %f %f %f %f', 'CommentStyle', '#', ...
        'Delimiter', {' ', '\t'}, 'MultipleDelimsAsOne', true);
    data.tow = c{1};
    data.dtheta = [c{2}, c{3}, c{4}];
    data.dvel = [c{5}, c{6}, c{7}];
    data.dt = [];
end

function data = readDecodedImuCsv(path)
    T = readtable(path);
    required = {'tow','dt','gyro_x','gyro_y','gyro_z','acc_x','acc_y','acc_z'};
    missing = required(~ismember(required, T.Properties.VariableNames));
    if ~isempty(missing)
        error('CSV missing columns: %s', strjoin(missing, ', '));
    end
    data.tow = T.tow;
    data.dt = T.dt;
    data.dtheta = [T.gyro_x, T.gyro_y, T.gyro_z];
    data.dvel = [T.acc_x, T.acc_y, T.acc_z];
end

function align = readCoarseAlignmentMeta(path)
    align.hasInfo = false;
    align.alignSeconds = NaN;
    align.imuStartTow = NaN;
    if ~exist(path, 'file')
        return;
    end
    fid = fopen(path, 'r');
    if fid < 0
        return;
    end
    cleanup = onCleanup(@() fclose(fid));
    while true
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
        line = strtrim(line);
        if startsWith(line, 'align_seconds')
            v = sscanf(line, 'align_seconds %f');
            if ~isempty(v), align.alignSeconds = v(1); end
        elseif startsWith(line, 'imu_start_tow')
            v = sscanf(line, 'imu_start_tow %f');
            if ~isempty(v), align.imuStartTow = v(1); end
        end
    end
    align.hasInfo = isfinite(align.alignSeconds) && isfinite(align.imuStartTow);
end

function shadePhase(tBoundary, colorVal, labelText)
    if ~isfinite(tBoundary) || tBoundary <= 0
        return;
    end
    ax = gca;
    yl = ylim(ax);
    p = patch(ax, [0 tBoundary tBoundary 0], [yl(1) yl(1) yl(2) yl(2)], colorVal, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.15);
    uistack(p, 'bottom');
    if ~isempty(labelText)
        text(ax, 0.985, 0.07, labelText, 'Units', 'normalized', ...
            'HorizontalAlignment', 'right', 'Color', [0.10 0.20 0.55], ...
            'FontSize', 9, 'FontWeight', 'bold');
    end
end
