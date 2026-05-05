clear; clc; close all;

%% Parameters
N_values = [20 40 60 80 100 150];
numWalks = 100;

mean_R = zeros(size(N_values));
std_R  = zeros(size(N_values));

mean_Rg = zeros(size(N_values));
std_Rg  = zeros(size(N_values));

trap_counts = zeros(size(N_values));
attempt_counts = zeros(size(N_values));
trap_rates = zeros(size(N_values));

blue = [0 0 1];

%% Run SAW simulations
for nIndex = 1:length(N_values)

    N = N_values(nIndex);

    R_values = zeros(numWalks,1);
    Rg_values = zeros(numWalks,1);

    successfulWalks = 0;
    trappedAttempts = 0;
    totalAttempts = 0;

    while successfulWalks < numWalks

        totalAttempts = totalAttempts + 1;

        [walk, trapped] = selfAvoidingWalk3D(N);

        if trapped
            trappedAttempts = trappedAttempts + 1;
        else
            successfulWalks = successfulWalks + 1;

            R_values(successfulWalks) = endToEndDistance(walk);
            Rg_values(successfulWalks) = radiusOfGyration(walk);
        end
    end

    mean_R(nIndex) = mean(R_values);
    std_R(nIndex)  = std(R_values);

    mean_Rg(nIndex) = mean(Rg_values);
    std_Rg(nIndex)  = std(Rg_values);

    trap_counts(nIndex) = trappedAttempts;
    attempt_counts(nIndex) = totalAttempts;
    trap_rates(nIndex) = trappedAttempts / totalAttempts;

    fprintf('N = %d completed\n', N);
    fprintf('  Successful walks = %d\n', successfulWalks);
    fprintf('  Trapped attempts = %d\n', trappedAttempts);
    fprintf('  Total attempts   = %d\n', totalAttempts);
    fprintf('  Trap rate        = %.3f\n\n', trap_rates(nIndex));
end

%% Example 3D self-avoiding walk
[exampleWalk, trapped] = selfAvoidingWalk3D(80);

while trapped
    [exampleWalk, trapped] = selfAvoidingWalk3D(80);
end

figure;
plot3(exampleWalk(:,1), exampleWalk(:,2), exampleWalk(:,3), '-o', ...
    'Color', blue, ...
    'MarkerFaceColor', blue, ...
    'LineWidth', 1.5);

grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
title('Example 3D Self-Avoiding Walk');

%% Bar graph: End-to-End Distance
figure;

b = bar(N_values, mean_R, 'FaceColor', blue);
b.FaceAlpha = 0.6;
hold on;

errorbar(N_values, mean_R, std_R, ...
    'k', ...
    'LineStyle', 'none', ...
    'LineWidth', 2, ...
    'CapSize', 12);

xlabel('Walk Length N');
ylabel('End-to-End Distance R');
title('SAW: Mean End-to-End Distance with Standard Deviation');
grid on;

%% Bar graph: Radius of Gyration
figure;

b = bar(N_values, mean_Rg, 'FaceColor', blue);
b.FaceAlpha = 0.6;
hold on;

errorbar(N_values, mean_Rg, std_Rg, ...
    'k', ...
    'LineStyle', 'none', ...
    'LineWidth', 2, ...
    'CapSize', 12);

xlabel('Walk Length N');
ylabel('Radius of Gyration R_g');
title('SAW: Mean Radius of Gyration with Standard Deviation');
grid on;

%% Trap counts
figure;

bar(N_values, trap_counts, 'FaceColor', [0.8 0.2 0.2]);

xlabel('Walk Length N');
ylabel('Number of Trapped Attempts');
title('Number of Trapped SAW Attempts');
grid on;

%% Trap rate
figure;

plot(N_values, trap_rates, '-o', ...
    'LineWidth', 2, ...
    'Color', [0.8 0.2 0.2], ...
    'MarkerFaceColor', [0.8 0.2 0.2]);

xlabel('Walk Length N');
ylabel('Trap Rate');
title('Trap Rate vs Walk Length');
grid on;

%% Log-log scaling plot
figure;

loglog(N_values, mean_R, '-o', ...
    'Color', [1,0,0], ...
    'LineWidth', 2);

hold on;

loglog(N_values, mean_Rg, '-s', ...
    'Color', blue, ...
    'LineWidth', 2);

xlabel('Walk Length N');
ylabel('Mean Size');
title('SAW Scaling Behavior');
legend('End-to-End Distance R', 'Radius of Gyration R_g');
grid on;

%% Scaling exponent
p_R  = polyfit(log(N_values), log(mean_R), 1);
p_Rg = polyfit(log(N_values), log(mean_Rg), 1);

fprintf('\nEstimated scaling exponent for R: %.3f\n', p_R(1));
fprintf('Estimated scaling exponent for Rg: %.3f\n', p_Rg(1));

fprintf('\nTrap statistics:\n');
for i = 1:length(N_values)
    fprintf('N = %d: trapped = %d, attempts = %d, trap rate = %.3f\n', ...
        N_values(i), trap_counts(i), attempt_counts(i), trap_rates(i));
end

%% ---------------- FUNCTIONS ----------------

function [walk, trapped] = selfAvoidingWalk3D(N)

    directions = [
         1  0  0;
        -1  0  0;
         0  1  0;
         0 -1  0;
         0  0  1;
         0  0 -1
    ];

    walk = zeros(N+1,3);
    walk(1,:) = [0 0 0];

    visited = containers.Map();
    visited('0_0_0') = true;

    trapped = false;

    for step = 1:N

        current = walk(step,:);
        order = randperm(6);
        moved = false;

        for j = 1:6

            dir = directions(order(j),:);
            new_pos = current + dir;

            key = sprintf('%d_%d_%d', ...
                new_pos(1), new_pos(2), new_pos(3));

            if ~isKey(visited, key)
                walk(step+1,:) = new_pos;
                visited(key) = true;
                moved = true;
                break;
            end
        end

        if ~moved
            trapped = true;
            walk = [];
            return;
        end
    end
end

function R = endToEndDistance(walk)
    R = norm(walk(end,:) - walk(1,:));
end

function Rg = radiusOfGyration(walk)
    r_cm = mean(walk,1);
    Rg = sqrt(mean(sum((walk - r_cm).^2,2)));
end
