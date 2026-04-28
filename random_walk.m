clear; clc; close all;


N = 100;              % number of monomers/steps
numWalks = 100;       % number of walks for Rg distribution

% Generate and visualize one self-avoiding walk
walk = selfAvoidingWalk3D(N);

figure;
plot3(walk(:,1), walk(:,2), walk(:,3), '-o', 'LineWidth', 1.5);
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
title('3D Self-Avoiding Random Walk: Uncharged Polymer');

%Rg plot
Rg_example = radiusOfGyration(walk);
fprintf('Example polymer Rg = %.4f\n', Rg_example);
Rg_values = zeros(numWalks,1);

for i = 1:numWalks
    walk_i = selfAvoidingWalk3D(N);
    Rg_values(i) = radiusOfGyration(walk_i);
end


figure;
histogram(Rg_values, 15);
xlabel('Radius of Gyration R_g');
ylabel('Frequency');
title('Distribution of Radius of Gyration for 100 Self-Avoiding Walks');
grid on;

% Summary statistics
fprintf('Mean Rg = %.4f\n', mean(Rg_values));
fprintf('Std Rg = %.4f\n', std(Rg_values));


%3D self-avoiding walk
function walk = selfAvoidingWalk3D(N)

    directions = [
        1  0  0;
       -1  0  0;
        0  1  0;
        0 -1  0;
        0  0  1;
        0  0 -1
    ];

    maxAttempts = 1000;

    for attempt = 1:maxAttempts

        walk = zeros(N,3);
        occupied = containers.Map;
        occupied("0,0,0") = true;

        success = true;

        for step = 2:N
            current = walk(step-1,:);
            order = randperm(6);
            moved = false;

            for j = 1:6
                newPos = current + directions(order(j),:);
                key = sprintf('%d,%d,%d', newPos(1), newPos(2), newPos(3));

                if ~isKey(occupied, key)
                    walk(step,:) = newPos;
                    occupied(key) = true;
                    moved = true;
                    break;
                end
            end

            if ~moved
                success = false;
                break;
            end
        end

        if success
            return;
        end
    end

    error('Failed to generate a self-avoiding walk. Try reducing N.');

end


% Function: Radius of Gyration
function Rg = radiusOfGyration(walk)

    centerOfMass = mean(walk, 1);

    squaredDistances = sum((walk - centerOfMass).^2, 2);

    Rg = sqrt(mean(squaredDistances));

end
