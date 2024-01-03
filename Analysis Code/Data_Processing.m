clear
clc

currentPath = pwd;
folders = split(currentPath, '\');
newPath = join(folders(1:length(folders)-2),"\");

addpath(strcat(newPath{1} , '\Stringer Data\Data'))
ExpData = matfile('natimg2800_M170717_MP033_2017-08-20.mat');

% Add Helper Function to Path
newPath = join(folders(1:length(folders)-1),"\");
addpath(strcat(newPath{1}, '\Helper Function'))



%%
Stim = ExpData.stim;
Response = Stim.resp;
Spont = Stim.spont;
Locations = ExpData.med;

offneuron = find(mean(Response)==0);
Response(:,offneuron) = [];
Locations(offneuron,:) = [];
Spont(:,offneuron) = [];

[normResponse, mu, sigma ]= zscore(Response,1,1);
[spontnormResponse, spontmu, spontsigma ]= zscore(Spont,1,1);
num_nuerons = size(Response,1);
%% View a random time series

randNueron = randi(num_nuerons);

plot(Response(1:30,randNueron))
hold on
plot((2*sigma(randNueron) + mu(randNueron)).*ones(30,1),'--', 'LineWidth',2)
title('Raw Resonse Time Series')
xlabel('Time')
ylabel('Response')
hold off


%% Histogram of data

histogram(normResponse, 'Normalization','pdf', 'FaceColor',[0,0,0])
hold on
for i = randi(num_nuerons, 1, 10)
    histogram(normResponse(:,i), 'FaceColor',[0.7,0.7,0.7], 'Normalization','pdf', 'FaceAlpha',0.3)
end
%set(gca ,'YScale', 'log')
xlim([-1,10])
xlabel('Z-score')
ylabel('Distribution')
legend('Aggregate Data', 'Random Neuron')
hold off


%% Binarizing Data
binary = (Response - mu)> 2*sigma;
%imshow(binary)


%% Filter Data
spiking_patterns = binary.';
num_bins = size(spiking_patterns,2);
num_nuerons = size(spiking_patterns,1);

datacorr = spiking_patterns*spiking_patterns.'/num_bins;
datamean = mean(spiking_patterns.');

datacorr_pseudo = datacorr + ones(num_nuerons)/(num_bins +1);

sum(datamean == 0)
Hind = sum(Entropy(datamean));

%% Minimal distance tree

[JTreedist, hTreedist, HTreedist] = Minimum_Distance_Tree(datamean,datacorr_pseudo,Locations);

%% Minimal distance GSP

[JGSPdist, hGSPdist, HGSPdist] = Minimum_Distance_GSP(datamean,datacorr_pseudo,Locations);

%% Find Optimal Tree

[JTree, hTree, HTree] = FindTree(datamean, datacorr_pseudo);

%% Random GSP

[JGSPrand, hGSPrand, HGSPrand] = RandomGSPFit(datamean, datacorr_pseudo);

%% Find Optimal GSP

% Ran on Cluster
tic
[hGSP, JGSP, HGSP, NewEnt2] = find_GSP_update(datamean, datacorr_pseudo);
toc

%% Entropy Drops

(1-HTreedist/Hind)*100
(1-HGSPdist/Hind)*100
(1-HGSPrand/Hind)*100
(1-HTree/Hind)*100
(1-HGSP/Hind)*100


%% Cumulative Spike Density 

numspikes = sum(full(binary));
[C,ia,ic] = unique(numspikes.','rows');
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];


figure
plot(C, a_counts/sum(a_counts), '.', 'MarkerSize', 10)

%% Plot Locations and network

xs = Locations(:,1);
ys = Locations(:,2);
zs = Locations(:,3);


zlayer = (zs >  0);

nodeLocations = Locations(zlayer,:);
adjacencyMatrix = (JGSPrand(zlayer,zlayer) ~= 0 );



% Plot nodes
scatter3(nodeLocations(:, 1), nodeLocations(:, 2), nodeLocations(:, 3), 'o', 'filled');
hold on;

% Plot edges
[numNodes, ~] = size(adjacencyMatrix);
for i = 1:numNodes
    for j = i+1:numNodes
        if adjacencyMatrix(i, j) == 1
            x = [nodeLocations(i, 1), nodeLocations(j, 1)];
            y = [nodeLocations(i, 2), nodeLocations(j, 2)];
            z = [nodeLocations(i, 3), nodeLocations(j, 3)];
            plot3(x, y, z, 'k-');
        end
    end
end

hold off;

% Set axis labels
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');

% Add a title
title('3D Graph Visualization');

% Adjust the view for better visibility
view(3);
grid on;


%% Distance

distanceMetric = squareform(pdist(Locations));
distanceMetric = distanceMetric(5006,:);

region   = (distanceMetric <  200);

reduceddatamean = datamean(region);
reduceddatacorr = datacorr_pseudo(region, region);
reducedLocation = Locations(region,:);


% Minimal distance tree

[JTreedist, hTreedist, HTreedist] = Minimum_Distance_Tree(reduceddatamean,reduceddatacorr,reducedLocation);

% Minimal distance GSP

[JGSPdist, hGSPdist, HGSPdist] = Minimum_Distance_GSP(reduceddatamean,reduceddatacorr,reducedLocation);

% Find Optimal Tree

[JTree, hTree, HTree] = FindTree(reduceddatamean, reduceddatacorr);

% Random GSP

[JGSPrand, hGSPrand, HGSPrand] = RandomGSPFit(reduceddatamean, reduceddatacorr);

% Find Optimal GSP

% Ran on Cluster
tic
[hGSP, JGSP, HGSP, NewEnt2] = find_GSP_update(reduceddatamean, reduceddatacorr);
toc

% Entropy Drops

(HTreedist/Hind)*100
(HGSPdist/Hind)*100
(HGSPrand/Hind)*100
(HTree/Hind)*100
(HGSP/Hind)*100

