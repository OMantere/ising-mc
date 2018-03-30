close all

results = load('results');
x = squeeze(results(:, 4));

figure
plot(x, results(:, 2))
title('Magnetization')
figure
plot(x, results(:, 3))
title('Heat capacity')
figure
plot(x, results(:, 1))
title('Energy')
