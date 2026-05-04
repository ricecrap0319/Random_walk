clear; clc; close all;

N = 100;

[walk, collisions, collisionCountPerStep] = randomWalk3D_withCollisionTracking(N);


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