clear; clc; close all;

N = 50;
nSteps = 50000;
sampleEvery = 100;

T = 1.0;             
epsilon = 0.5;        
sigma = 1.0;
k_spring = 50;
r0 = 1.0;

moveSize = 0.25;

chain = zeros(N,3);
for i = 1:N
    chain(i,:) = [(i-1)*r0, 0, 0];
end


E_current = polymerEnergy(chain, epsilon, sigma, k_spring, r0);

Rg_values = [];
E_values = [];

accepted = 0;
totalMoves = 0;


for step = 1:nSteps

    old_chain = chain;
    old_E = E_current;

    idx = randi(N);

    displacement = moveSize * randn(1,3);
    chain(idx,:) = chain(idx,:) + displacement;
    E_new = polymerEnergy(chain, epsilon, sigma, k_spring, r0);
    dE = E_new - old_E;

    totalMoves = totalMoves + 1;

    if dE <= 0 || rand < exp(-dE / T)
        E_current = E_new;
        accepted = accepted + 1;
    else
        chain = old_chain;
        E_current = old_E;
    end
    if mod(step, sampleEvery) == 0
        Rg_values(end+1) = radiusOfGyration(chain);
        E_values(end+1) = E_current;
    end
end

acceptanceRate = accepted / totalMoves;
fprintf('Acceptance rate = %.3f\n', acceptanceRate);

figure;
plot3(chain(:,1), chain(:,2), chain(:,3), '-o','LineWidth',2);
grid on; axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
title('Final Polymer Configuration');


figure;
histogram(Rg_values, 25, 'Normalization','pdf');
xlabel('R_g'); ylabel('P(R_g)');
title('Radius of Gyration Distribution');


[counts, edges] = histcounts(Rg_values,25,'Normalization','pdf');
Rg_centers = (edges(1:end-1)+edges(2:end))/2;

F = -T * log(counts + 1e-12);
F = F - min(F);

figure;
plot(Rg_centers, F, '-o','LineWidth',2);
xlabel('R_g'); ylabel('F(R_g)');
title('Free Energy Profile');
grid on;


figure;
plot(E_values,'LineWidth',2);
xlabel('Sample'); ylabel('Energy');
title('Energy Evolution');
grid on;


function Rg = radiusOfGyration(chain)
    r_cm = mean(chain,1);
    Rg = sqrt(mean(sum((chain - r_cm).^2,2)));
end

function E = polymerEnergy(chain, epsilon, sigma, k_spring, r0)

    N = size(chain,1);
    E = 0;

   
    for i = 1:N-1
        r = norm(chain(i+1,:) - chain(i,:));
        E = E + 0.5 * k_spring * (r - r0)^2;
    end

    for i = 1:N
        for j = i+2:N

            r = norm(chain(i,:) - chain(j,:));

            if r < 1e-6
                E = E + 1e12;
            else
                E = E + 4 * epsilon * ((sigma/r)^12 - (sigma/r)^6);
            end

        end
    end
end
