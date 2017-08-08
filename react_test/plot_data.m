clear variables
close all

concfile = fopen('time_concs.txt', 'r');
shape_concs = fscanf(concfile,'%f', 3);
concs = fscanf(concfile,'%f');
concs = reshape(concs, shape_concs');
fclose(concfile);

x = linspace(0, 0.5, shape_concs(1));

for i = 1 : shape_concs(3)

    figure(1)

    plot(x, concs(:, 2, i))
%     title(['t=',num2str(i),'*dt'])
%     set(gca,'ylim',[0,0.36867]);

    pause(0.02)
end