clear; clc; close all;

N = 100;
numWalks = 100;

[walk, collisions, collisionCountPerStep] = randomWalk3D_withCollisionTracking(N);
Rg_example = radiusOfGyration(walk);


figure;

subplot(1,2,1);

plot3(walk(:,1), walk(:,2), walk(:,3), '-o', 'LineWidth', 1.5);
grid on;
axis equal;
xlabel('x'); ylabel('y'); zlabel('z');

title(sprintf('3D Random Walk (%d collisions)', collisions));


subplot(1,2,2);

plot(1:N, collisionCountPerStep, 'LineWidth', 2);
xlabel('Step');
ylabel('Cumulative Collisions');

title('Collision Growth Over Time');
grid on;

fprintf('Example polymer Rg = %.4f\n', Rg_example);

Rg_values = zeros(numWalks,1);

for i = 1:numWalks
    walk_i = randomWalk3D(N);
    Rg_values(i) = radiusOfGyration(walk_i);
end

figure;
histogram(Rg_values, 15);
xlabel('Radius of Gyration R_g');
ylabel('Frequency');
title('Distribution of Radius of Gyration for 100 Non-Self-Avoiding Walks');
grid on;

fprintf('Mean Rg = %.4f\n', mean(Rg_values));
fprintf('Std Rg = %.4f\n', std(Rg_values));

function walk = randomWalk3D(N)

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

    for step = 1:N
        dir = directions(randi(6),:);
        walk(step+1,:) = walk(step,:) + dir;
    end
end

function [walk, collisions, collisionCountPerStep] = randomWalk3D_withCollisionTracking(N)

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

    collisions = 0;
    collisionCountPerStep = zeros(N,1);

    for step = 1:N
        dir = directions(randi(6),:);
        new_pos = walk(step,:) + dir;

        key = sprintf('%d_%d_%d', new_pos(1), new_pos(2), new_pos(3));

        if isKey(visited, key)
            collisions = collisions + 1;
        else
            visited(key) = true;
        end

        collisionCountPerStep(step) = collisions;

        walk(step+1,:) = new_pos;
    end
end

function Rg = radiusOfGyration(walk)

    centerOfMass = mean(walk, 1);
    squaredDistances = sum((walk - centerOfMass).^2, 2);
    Rg = sqrt(mean(squaredDistances));

end
