
Cost_G = xlsread('Cost_G.csv');
Diam_G = xlsread('Diam_G.csv');
Tanks_G = xlsread('Tanks_G.csv');


figure('Color','w')
image(Cost_G,'CDataMapping','scaled')
xlabel('Inital tank radius (0.12-2.5) [m]')
ylabel('Inital number of tanks (1-25) [m]')
title(' SENSITIVITY : COST OF ENERY THERMAL STORAGE')
colormap('jet')
axis tight

figure('Color','w')
image(Diam_G,'CDataMapping','scaled')
xlabel('Inital tank radius (0.12-2.5) [m]')
ylabel('Inital number of tanks (1-25) [m]')
title(' SENSITIVITY : TANKS DIAMETER')
colormap('jet')
axis tight

figure('Color','w')
image(Tanks_G,'CDataMapping','scaled')
xlabel('Inital tank radius (0.12-2.5) [m]')
ylabel('Inital number of tanks (1-25) [m]')
title(' SENSITIVITY : NUMBER OF TANKS')
colormap('jet')
axis tight