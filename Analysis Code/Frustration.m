%% Frustration Analysis
clear
clc

currentPath = pwd;
folders = split(currentPath, '\');
newPath = join(folders(1:length(folders)-2),"\");

addpath(strcat(newPath{1} , '\Model Data'))

load("natimg2800_2017-08-20_Fitting.mat")


% Add Helper Function to Path
newPath = join(folders(1:length(folders)-1),"\");
addpath(strcat(newPath{1}, '\Helper Function'))


%% Optimal

[G,D] = decimate_GSP(JGSP);

frustrations = zeros(num_nuerons-2,1);

count = 1;
for keep = D(1:size(D,1)-1,:).'
    
    
    frustrations(count) = JGSP(keep(1),keep(2))*JGSP(keep(2),keep(3))*JGSP(keep(3),keep(1));

    count = count + 1;

    
end

histogram(frustrations, 'Normalization','probability')
hold on

% Random

[G,D] = decimate_GSP(JGSPrand);

frustrations = zeros(num_nuerons-2,1);

count = 1;
for keep = D(1:size(D,1)-1,:).'
    
    
    frustrations(count) = JGSPrand(keep(1),keep(2))*JGSPrand(keep(2),keep(3))*JGSPrand(keep(3),keep(1));

    count = count + 1;

    
end

histogram(frustrations, 'Normalization','probability')

% Min Dist

[G,D] = decimate_GSP(JGSPdist);

frustrations = zeros(num_nuerons-2,1);

count = 1;
for keep = D(1:size(D,1)-1,:).'
    
    
    frustrations(count) = JGSPdist(keep(1),keep(2))*JGSPdist(keep(2),keep(3))*JGSPdist(keep(3),keep(1));

    count = count + 1;

    
end

histogram(frustrations, 'Normalization','probability')

hold off
legend('Optimal','Random','Min Dist')
xlim([-15,15])
