clear; clc; close all;

%% Parameters
N_values = [20 40 60 80 100 150 200];
numWalks = 200;

mean_R = zeros(size(N_values));
std_R  = zeros(size(N_values));

mean_Rg = zeros(size(N_values));
std_Rg  = zeros(size(N_values));

%% Simulations
for nIndex = 1:length(N_values)

    N = N_values(nIndex);

    R_values = zeros(numWalks,1);
    Rg_values = zeros(numWalks,1);

    for i = 1:numWalks
        walk = randomWalk3D(N);

        R_values(i) = endToEndDistance(walk);
        Rg_values(i) = radiusOfGyration(walk);
    end

    mean_R(nIndex) = mean(R_values);
    std_R(nIndex)  = std(R_values);

    mean_Rg(nIndex) = mean(Rg_values);
    std_Rg(nIndex)  = std(Rg_values);

    fprintf('N = %d completed\n', N);
end

%% Example 3D walk
exampleWalk = randomWalk3D(100);

figure;
plot3(exampleWalk(:,1), exampleWalk(:,2), exampleWalk(:,3), '-o', 'LineWidth', 1.5);
grid on;
axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
title('Example 3D Non-Self-Avoiding Random Walk');

%% -------- BAR GRAPH: End-to-End Distance --------
figure;
bar(N_values, mean_R);
hold on;
errorbar(N_values, mean_R, std_R, 'k.', 'LineWidth', 1.5);

xlabel('Walk Length N');
ylabel('End-to-End Distance R');
title('Mean End-to-End Distance (Random Walk)');
grid on;

%% -------- BAR GRAPH: Radius of Gyration --------
figure;
bar(N_values, mean_Rg);
hold on;
errorbar(N_values, mean_Rg, std_Rg, 'k.', 'LineWidth', 1.5);

xlabel('Walk Length N');
ylabel('Radius of Gyration R_g');
title('Mean Radius of Gyration (Random Walk)');
grid on;

%% -------- LOG-LOG SCALING --------
figure;
loglog(N_values, mean_R, '-o', 'LineWidth', 2);
hold on;
loglog(N_values, mean_Rg, '-s', 'LineWidth', 2);

xlabel('Walk Length N');
ylabel('Mean Size');
title('Scaling Behavior (Random Walk)');
legend('R', 'R_g');
grid on;

%% Estimate scaling exponent
p_R  = polyfit(log(N_values), log(mean_R), 1);
p_Rg = polyfit(log(N_values), log(mean_Rg), 1);

fprintf('\nScaling exponent (R):  %.3f\n', p_R(1));
fprintf('Scaling exponent (Rg): %.3f\n', p_Rg(1));

%% -------- FUNCTIONS --------

function walk = randomWalk3D(N)

    directions = [
         1  0  0;
        -1  0  0;
         0  1  0;
         0 -1  0;
         0  0  1;
         0  0 -1
    ];

    walk = zeros(N+1, 3);
    walk(1,:) = [0 0 0];

    for step = 1:N
        dir = directions(randi(6), :);
        walk(step+1,:) = walk(step,:) + dir;
    end
end

function R = endToEndDistance(walk)
    R = norm(walk(end,:) - walk(1,:));
end

function Rg = radiusOfGyration(walk)
    r_cm = mean(walk, 1);
    Rg = sqrt(mean(sum((walk - r_cm).^2, 2)));
end
