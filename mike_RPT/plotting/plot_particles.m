clear variables
close all

concfile = fopen('time_concs.txt', 'r');
num = fscanf(concfile,'%f', 1);


% 
% x = linspace(0, 0.5, shape_concs(1));
% 
for i = 1 : num

    figure(1)
    
    nactive = fscanf(concfile,'%f', 1);
    locs = fscanf(concfile,'%f', nactive);
    concA = fscanf(concfile,'%f', nactive);
    concB = fscanf(concfile,'%f', nactive);

    scatter(locs, concA)
    hold on
    scatter(locs, concB)
%     title(['t=',num2str(i),'*dt'])
	axis([0, 0.5, 0, 2])
    hold off

    pause(0.02)
end

fclose(concfile);

figure(2)
hist(locs,100)