%%
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


datacorr = datacorr_pseudo;

hind = log(datamean);
[B,I] = sort(datamean);

binary_M = binary.';
%%


simul = zeros(num_nuerons,num_nuerons);

for t = 1:size(binary_M,2)
    fire = find(binary_M(:,t));
    simul(fire,fire) = ones(length(fire));
end


%% Plot 2 d
figure
plot(datamean(I), hind(I), ':', LineWidth=2, Color='black')
hold on
plot(datamean,hTree,'.', MarkerSize=14)




plot(datamean,hGSP,'.', MarkerSize=14, Color='blue')
plot(datamean,hGSPdist,'.', MarkerSize=14, Color=[0.2,0.8,0])
legend('Independent', 'Tree', 'GSP','Min Dist.', 'Location','northwest','FontSize',14)
xlabel('Average activity $\left<x_i\right>$', 'Interpreter','latex', 'FontSize',14)
ylabel('External field $h_i$', 'Interpreter','latex','FontSize',14)
hold off

%% Plot 2 e

Jtemp = JTree;
Jtemp(simul == 0) = 0;
figure
plot(hTree, sum(Jtemp),'.', MarkerSize=10)
percentTree = sum(sum(Jtemp) >= -hTree)/num_nuerons
hold on
Jtemp = JGSP;
Jtemp(simul == 0) = 0;
plot(hGSP,sum(Jtemp),'.', MarkerSize=10)
percentGSP = sum(sum(Jtemp) >= -hGSP)/num_nuerons

Jtemp = JGSPdist;
Jtemp(simul == 0) = 0;
plot(hGSP,sum(Jtemp),'.', MarkerSize=10)
percentGSPdist = sum(sum(Jtemp) >= -hGSP)/num_nuerons
plot(linspace(-12,-1),-linspace(-12,-1),':',LineWidth=2, Color='black')

plot(linspace(-12,-1),-linspace(-12,-1),':',LineWidth=2, Color='black')
ylabel('Total interaction $\sum_j J_{ij}$', 'Interpreter','latex','FontSize',14)
xlabel('External field $h_i$', 'Interpreter','latex','FontSize',14)
legend(sprintf('Tree %.2f %% ', percentTree*100 ), sprintf('GSP %.2f %%', percentGSP*100 ),sprintf('GSP min dist. %.2f %%', percentGSPdist*100 ), 'y = -x','FontSize',14)



%% Plot 2 b 

figure
t = tiledlayout(1,1);

ax1 = axes(t);
barh([(Hind - HGSP)/Hind,(Hind - HTree)/Hind, (Hind - HGSPdist)/Hind ])
xlabel(ax1, "I^T/S^{ind}")
xlim(ax1,[0,200/Hind])
yticklabels({'GSP', 'Tree', 'Min Dist. GSP'})
ax2 = axes(t);
barh([ (Hind - HGSP)/log(2), (Hind - HTree)/log(2),(Hind - HGSPdist)/log(2)])
ax2.XAxisLocation = 'top';
xlabel(ax2, "I^T (bits)")
xlim(ax2,[0,200/log(2)])
yticklabels({' ',' ', ' '})



%% Plot 2 c

figure
Jtemp = JTree;
Jtemp(simul == 0) = 0;
histogram(Jtemp(Jtemp ~= 0),20,'Normalization','pdf');
hold on
Jtemp = JGSP;
Jtemp(simul == 0) = 0;
histogram(Jtemp(Jtemp ~= 0),20, 'Normalization','pdf');
hold off
ylabel('Prob. density', 'FontSize',14)
xlabel('Interactions J_{ij}', 'FontSize',14)
legend('Tree','GSP','FontSize',14)


%% Plot 3 Synchrony

clc
numspikes = sum(binary_M);
[C,ia,ic] = unique(numspikes.','rows');
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];
sum(numspikes)



[~, randPerms] = sort(rand(size(binary_M)));


shuffled = zeros(size(binary_M));

for neuron = 1:num_nuerons

    shuffled(neuron,:) = binary_M(neuron, randperm(num_bins));
end

numspikes = sum(shuffled);
[Cshuffled,ia,ic] = unique(numspikes.','rows');
a_countsshuffled = accumarray(ic,1);
value_counts = [Cshuffled, a_countsshuffled];



sum(numspikes)


[spikesGSP, probspikeGSP] = Estimate_NumSpikeDist(JGSP,hGSP,1);
[spikesTree, probspikeTree] = Estimate_NumSpikeDist(JTree,hTree,1);

%% continue Plot 3

figure
semilogy(C, a_counts/sum(a_counts), '.', 'MarkerSize', 10)
hold on
semilogy(Cshuffled, a_countsshuffled/sum(a_countsshuffled), '.', 'MarkerSize', 10)

semilogy(spikesGSP,probspikeGSP, '.', 'MarkerSize', 10)
semilogy(spikesTree,probspikeTree, '.', 'MarkerSize', 10)
title('Spiking Probablity')
ylabel('Probability')
xlabel('Number of Spikes')
legend('Observed','Shuffled','GSP', 'Tree')
xlim([0,1000])
hold off

%% Scatter Plot of Mutual Information


Mutual_info = MI2(datamean,datacorr);
psueddatamean = datamean + 1/num_bins;

datacorr = datacorr_pseudo - eye(num_nuerons)/(num_bins +1);

corr_coef = (datacorr - datamean.*datamean.')./sqrt((datamean-datamean.^2).*(datamean-datamean.^2).');

%corr_coef(simul == 0) = nan;
%%
figure
scatter(reshape(corr_coef,[1,num_nuerons*num_nuerons]), reshape(Mutual_info,[1,num_nuerons*num_nuerons]), '.')
xlim([-0.15,0.85]) 
xlabel('Correlation coefficient')
ylabel('Mutual info (bits)')

%% Histogram of Correlation Coeficients

figure
histogram(corr_coef(~isnan(corr_coef)),40,'Normalization','pdf');
set(gca,'YScale','log')
xlim([-0.15,0.85])
ylim([10^(-4),10^2])
xlabel('Correlation coefficient','FontSize',14)
ylabel('Prob. density','FontSize',14)

%% Probability Distribution of Mutual Information

MIS2 = reshape(Mutual_info,[1,num_nuerons*num_nuerons]);
MIS2 = sort(MIS2);
MIS2(isnan(MIS2)) = 0;

numbins = 1000;
[Y, E] = discretize(MIS2,numbins);

[C,ia,ic] = unique(Y);
a_counts = accumarray(ic,1);

figure
loglog(E(Y(ia))/log(2), a_counts/sum(a_counts), '.', MarkerSize=10)
xlabel('Mutual Information I_{ij} (bits)')
ylabel('Probability')


%% Histogram of Topologies


distanceMetric = squareform(pdist(Locations));
bins = logspace(log10(1),log10(1500),100);

figure;
bJ = JGSP ~= 0;
bJ = double(bJ);
bJ(bJ == 0) = nan;
histogram(distanceMetric.*bJ, bins,Normalization="probability", FaceColor=[1,0,0])
hold on

bins = logspace(log10(1),log10(1500),100);
bJ = JGSPrand ~= 0;
bJ = double(bJ);
bJ(bJ == 0) = nan;
histogram(distanceMetric.*bJ,bins, Normalization="probability",FaceAlpha=0.3, FaceColor=[0,1,0])

bins = logspace(log10(1),log10(1500),100);
bJ = JGSPdist ~= 0;
bJ = double(bJ);
bJ(bJ == 0) = nan;
histogram(distanceMetric.*bJ, bins, Normalization="probability",FaceAlpha=0.3, FaceColor=[0,0,1])

bins = logspace(log10(1),log10(1500),100);
bJ = JTree ~= 0;
bJ = double(bJ);
bJ(bJ == 0) = nan;
histogram(distanceMetric.*bJ,bins, Normalization="probability")

hold off

legend('Optimal GSP', 'GSP Rand', 'GSP Min Dist.', 'Optimal Tree')
xlabel('Pairwise Distance')
ylabel('Probability')
set(gca, "XScale", "log")

%%

[modelmean, modelcorr] = findStats(JGSP,hGSP);

model_coef = (modelcorr - modelmean.*modelmean.')./sqrt((modelmean-modelmean.^2).*(modelmean-modelmean.^2).');

figure
plot(corr_coef(JGSP ~= 0),model_coef(JGSP ~= 0), '.')
hold on
plot(linspace(-0.1,0.9),linspace(-0.1,0.9),'--', Color='black')
hold off
xlabel('Data Corr Coef', 'FontSize',16)
ylabel('Model Corr Coef', 'FontSize',16)

%%

[m, C, X, Z] = correlations_GSP_01(JGSP, hGSP);
%%


model_coef = X./sqrt((m-m.^2).*(m-m.^2).');

figure
plot(corr_coef(JGSP ~= 0 ),model_coef(JGSP ~= 0), '.')
hold on
plot(corr_coef(metrics(:,:,10) ==  2),model_coef(metrics(:,:,10) == 2), '.')
plot(corr_coef(metrics(:,:,10) ==  3),model_coef(metrics(:,:,10) == 3), '.')
plot(linspace(-0.1,0.9),linspace(-0.1,0.9),'--', Color='black')
hold off
xlabel('Data Corr Coef', 'FontSize',16)
ylabel('Model Corr Coef', 'FontSize',16)
legend('Top. dis. 1','Top. dis. 2','Top. dis. 3','Location','northwest', 'FontSize',14)

%% Plots for 1,2,3 topological distance


plotvariance(corr_coef(JGSP ~= 0 ),model_coef(JGSP ~= 0), 100, 'blue')
hold on
plotvariance(corr_coef(triu(metrics(:,:,10) ==  2)), model_coef(triu(metrics(:,:,10) == 2)), 11 , 'red')
plotvariance(corr_coef(triu(metrics(:,:,10) ==  3)), model_coef(triu(metrics(:,:,10) == 3)), 6, [0.9290 0.6940 0.1250])
plot(linspace(-0.1,0.9),linspace(-0.1,0.9),'--', Color='black')
hold off
xlabel('Data Corr Coef', 'FontSize',16)
ylabel('Model Corr Coef', 'FontSize',16)
legend('','Top. dis. 1','','Top. dis. 2','','Top. dis. 3','Location','northwest', 'FontSize',14)

%%
plot(linspace(-0.1,0.9),linspace(-0.1,0.9),'--', Color='black')
hold on
plotvariance(corr_coef(triu(eye(num_nuerons)~=1)), model_coef(triu(eye(num_nuerons)~=1)), 11 , 'red')
xlabel('Data Corr Coef', 'FontSize',16)
ylabel('Model Corr Coef', 'FontSize',16)
title('Full Data')

%% Optimal Tree

[m, C, X, Z] = correlations_GSP_01(JTree, hTree);
%


model_coef = X./sqrt((m-m.^2).*(m-m.^2).');

figure
plot(corr_coef(JTree ~= 0 ),model_coef(JTree ~= 0), '.')
hold on
plot(corr_coef(metrics(:,:,8) ==  2),model_coef(metrics(:,:,8) == 2), '.')
plot(corr_coef(metrics(:,:,8) ==  3),model_coef(metrics(:,:,8) == 3), '.')
plot(linspace(-0.1,0.9),linspace(-0.1,0.9),'--', Color='black')
hold off
xlabel('Data Corr Coef', 'FontSize',16)
ylabel('Model Corr Coef', 'FontSize',16)
legend('Top. dis. 1','Top. dis. 2','Top. dis. 3','Location','northwest', 'FontSize',14)

% Plots for 1,2,3 topological distance

figure
plotvariance(corr_coef(JTree ~= 0 ),model_coef(JTree ~= 0), 100, 'blue')
hold on
plotvariance(corr_coef(triu(metrics(:,:,8) ==  2)), model_coef(triu(metrics(:,:,8) == 2)), 11 , 'red')
plotvariance(corr_coef(triu(metrics(:,:,8) ==  3)), model_coef(triu(metrics(:,:,8) == 3)), 6, [0.9290 0.6940 0.1250])
plot(linspace(-0.1,0.9),linspace(-0.1,0.9),'--', Color='black')
hold off
xlabel('Data Corr Coef', 'FontSize',16)
ylabel('Model Corr Coef', 'FontSize',16)
legend('','Top. dis. 1','','Top. dis. 2','','Top. dis. 3','Location','northwest', 'FontSize',14)
%%
figure 
plot(linspace(-0.1,0.9),linspace(-0.1,0.9),'--', Color='black')
hold on
plotvariance(corr_coef(triu(eye(num_nuerons)~=1)), model_coef(triu(eye(num_nuerons)~=1)), 11 , 'red')
xlabel('Data Corr Coef', 'FontSize',16)
ylabel('Model Corr Coef', 'FontSize',16)
title('Full Data Optimal Tree')


%% Min Dist GSP

[m, C, X, Z] = correlations_GSP_01(JGSPdist, hGSPdist);
%


model_coef = X./sqrt((m-m.^2).*(m-m.^2).');

figure
plot(corr_coef(JGSPdist ~= 0 ),model_coef(JGSPdist ~= 0), '.')
hold on
plot(corr_coef(metrics(:,:,6) ==  2),model_coef(metrics(:,:,6) == 2), '.')
plot(corr_coef(metrics(:,:,6) ==  3),model_coef(metrics(:,:,6) == 3), '.')
plot(linspace(-0.1,0.9),linspace(-0.1,0.9),'--', Color='black')
hold off
xlabel('Data Corr Coef', 'FontSize',16)
ylabel('Model Corr Coef', 'FontSize',16)
legend('Top. dis. 1','Top. dis. 2','Top. dis. 3','Location','northwest', 'FontSize',14)
title('Min Dist. GSP')

% Plots for 1,2,3 topological distance

figure
plotvariance(corr_coef(JGSPdist ~= 0 ),model_coef(JGSPdist ~= 0), 100, 'blue')
hold on
plotvariance(corr_coef(triu(metrics(:,:,6) ==  2)), model_coef(triu(metrics(:,:,6) == 2)), 11 , 'red')
plotvariance(corr_coef(triu(metrics(:,:,6) ==  3)), model_coef(triu(metrics(:,:,6) == 3)), 6, [0.9290 0.6940 0.1250])
plot(linspace(-0.1,0.9),linspace(-0.1,0.9),'--', Color='black')
hold off
xlabel('Data Corr Coef', 'FontSize',16)
ylabel('Model Corr Coef', 'FontSize',16)
legend('','Top. dis. 1','','Top. dis. 2','','Top. dis. 3','Location','northwest', 'FontSize',14)
title('Min Dist. GSP')
%%
figure 
plot(linspace(-0.1,0.9),linspace(-0.1,0.9),'--', Color='black')
hold on
plotvariance(corr_coef(triu(eye(num_nuerons)~=1)), model_coef(triu(eye(num_nuerons)~=1)), 11 , 'red')
xlabel('Data Corr Coef', 'FontSize',16)
ylabel('Model Corr Coef', 'FontSize',16)
title('Full Data Optimal Min Dist GSP')

%% Min Dist Tree

[m, C, X, Z] = correlations_GSP_01(JTreedist, hTreedist);
%


model_coef = X./sqrt((m-m.^2).*(m-m.^2).');

figure
plot(corr_coef(JTreedist ~= 0 ),model_coef(JTreedist ~= 0), '.')
hold on
plot(corr_coef(metrics(:,:,4) ==  2),model_coef(metrics(:,:,4) == 2), '.')
plot(corr_coef(metrics(:,:,4) ==  3),model_coef(metrics(:,:,4) == 3), '.')
plot(linspace(-0.1,0.9),linspace(-0.1,0.9),'--', Color='black')
hold off
xlabel('Data Corr Coef', 'FontSize',16)
ylabel('Model Corr Coef', 'FontSize',16)
legend('Top. dis. 1','Top. dis. 2','Top. dis. 3','Location','northwest', 'FontSize',14)
title('Min Dist. Tree')

% Plots for 1,2,3 topological distance

figure
plotvariance(corr_coef(JTreedist ~= 0 ),model_coef(JTreedist ~= 0), 100, 'blue')
hold on
plotvariance(corr_coef(triu(metrics(:,:,4) ==  2)), model_coef(triu(metrics(:,:,4) == 2)), 11 , 'red')
plotvariance(corr_coef(triu(metrics(:,:,4) ==  3)), model_coef(triu(metrics(:,:,4) == 3)), 6, [0.9290 0.6940 0.1250])
plot(linspace(-0.1,0.9),linspace(-0.1,0.9),'--', Color='black')
hold off
xlabel('Data Corr Coef', 'FontSize',16)
ylabel('Model Corr Coef', 'FontSize',16)
legend('','Top. dis. 1','','Top. dis. 2','','Top. dis. 3','Location','northwest', 'FontSize',14)
title('Min Dist. Tree')
%%
figure 
plot(linspace(-0.1,0.9),linspace(-0.1,0.9),'--', Color='black')
hold on
plotvariance(corr_coef(triu(eye(num_nuerons)~=1)), model_coef(triu(eye(num_nuerons)~=1)), 11 , 'red')
xlabel('Data Corr Coef', 'FontSize',16)
ylabel('Model Corr Coef', 'FontSize',16)
title('Full Data Optimal Min Dist Tree')

%%

%% GSP rand

[m, C, X, Z] = correlations_GSP_01(JGSPrand, hGSPrand);
%


model_coef = X./sqrt((m-m.^2).*(m-m.^2).');

figure
plot(corr_coef(JGSPrand ~= 0 ),model_coef(JGSPrand ~= 0), '.')
hold on
plot(corr_coef(metrics(:,:,2) ==  2),model_coef(metrics(:,:,2) == 2), '.')
plot(corr_coef(metrics(:,:,2) ==  3),model_coef(metrics(:,:,2) == 3), '.')
plot(linspace(-0.1,0.9),linspace(-0.1,0.9),'--', Color='black')
hold off
xlabel('Data Corr Coef', 'FontSize',16)
ylabel('Model Corr Coef', 'FontSize',16)
legend('Top. dis. 1','Top. dis. 2','Top. dis. 3','Location','northwest', 'FontSize',14)
title('Random GSP')

% Plots for 1,2,3 topological distance

figure
plotvariance(corr_coef(JGSPrand ~= 0 ),model_coef(JGSPrand   ~= 0), 100, 'blue')
hold on
plotvariance(corr_coef(triu(metrics(:,:,2) ==  2)), model_coef(triu(metrics(:,:,2) == 2)), 11 , 'red')
plotvariance(corr_coef(triu(metrics(:,:,2) ==  3)), model_coef(triu(metrics(:,:,2) == 3)), 6, [0.9290 0.6940 0.1250])
plot(linspace(-0.1,0.9),linspace(-0.1,0.9),'--', Color='black')
hold off
xlabel('Data Corr Coef', 'FontSize',16)
ylabel('Model Corr Coef', 'FontSize',16)
legend('','Top. dis. 1','','Top. dis. 2','','Top. dis. 3','Location','northwest', 'FontSize',14)
title('Random GSP')
%%
figure 
plot(linspace(-0.1,0.9),linspace(-0.1,0.9),'--', Color='black')
hold on
plotvariance(corr_coef(triu(eye(num_nuerons)~=1)), model_coef(triu(eye(num_nuerons)~=1)), 11 , 'red')
xlabel('Data Corr Coef', 'FontSize',16)
ylabel('Model Corr Coef', 'FontSize',16)
title('Full Data Random GSP')

%% Triplets

tic
[G,D] = decimate_GSP(JGSP);
toc

Data = zeros(100,2);
count = 1;
for keep = D(1:2,:).'

    tic
    [Gtest, Dtest, Jtest, htest] = decimate_GSP_carefully(JGSP,hGSP, keep);
    toc

    [meanspin, corr,threepoint] = Exact_Ising(Jtest(keep,keep), htest(keep).', 1);
    
    Data(1,count) = threepoint(1,2,3);
    
    Data(2,count) = (binary_M(keep(1),:).*binary_M(keep(2),:))*binary_M(keep(3),:).'/size(binary_M,2);
    count = count + 1

    Data(1,count)
    Data(2,count)
end

%%
plot(Data(:,1),Data(:,2), '.')
hold on
plot([min(Data(:,1)),max(Data(:,1))], [min(Data(:,1)),max(Data(:,1))])
hold off
xlabel('<x_ix_jx_k> Data')
ylabel('<x_ix_jx_k> Model')

%%

plotvariance(Data(:,1),Data(:,2), 10, 'blue')
hold on
plot([min(Data(:,1)),max(Data(:,1))], [min(Data(:,1)),max(Data(:,1))])
hold off
xlabel('<x_ix_jx_k> Data')
ylabel('<x_ix_jx_k> Model')

%% Cumulants

count = 1;
cumData = zeros(num_nuerons ,2);
for keep = D(:,:).'

    cumData(count,1) = Data(count, 1) - datamean(keep(1))*datacorr(keep(2),keep(3)) - datamean(keep(2))*datacorr(keep(3),keep(1))-datamean(keep(3))*datacorr(keep(1),keep(2)) + 2*datamean(keep(1))*datamean(keep(2))*datamean(keep(3));
    cumData(count,2) = Data(count, 2) - datamean(keep(1))*datacorr(keep(2),keep(3)) - datamean(keep(2))*datacorr(keep(3),keep(1))-datamean(keep(3))*datacorr(keep(1),keep(2)) + 2*datamean(keep(1))*datamean(keep(2))*datamean(keep(3));
    
    count = count + 1
 
end
%%

plot(cumData(:,1),cumData(:,2), '.')
hold on
plot([min(cumData(:,1)),max(cumData(:,1))], [min(cumData(:,1)),max(cumData(:,1))])
hold off
xlabel('[x_ix_jx_k] Data')
ylabel('[x_ix_jx_k] Model')

%%

plotvariance(cumData(:,1),cumData(:,2), 10, 'blue')
hold on
plot([min(cumData(:,1)),max(cumData(:,1))], [min(cumData(:,1)),max(cumData(:,1))])
hold off
xlabel('[x_ix_jx_k] Data')
ylabel('[x_ix_jx_k] Model')


%% Distance

distanceMetric = squareform(pdist(Locations));
distanceMetric = distanceMetric(5006,:);

region   = (distanceMetric <  30);

reduceddatamean = datamean(region);
reduceddatacorr = datacorr_pseudo(region);




%% Plot Locations and network

xs = Locations(:,1);
ys = Locations(:,2);
zs = Locations(:,3);

distanceMetric = squareform(pdist(Locations));
distanceMetric = distanceMetric(5006,:);



zlayer = (distanceMetric <  30);

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


% function [G, D, Jeff, heff] = decimate_GSP_carefully(G0, h, keep)
%     Input: nxn unweighted adjacency matrix G0 representing a GSP network.
% 
%     Output: Randomly remove nodes of degree 1 or 2. If we remove a node of
%     degree 2 then we place a new edge between the two neighbors of the
%     removed node. Continue this process until no more nodes of degree 1 or 2
%     exist and return the final adjacency matrix G. Also return the order of
%     decimations in the numDec x 3 matrix D, where D(t,1) is the node removed
%     at step t and D(t,2) and D(t,3) are the two nodes from which the node was
%     removed. We note that D(t,3) = 0 if the removed node had degree 1.
% 
%     Size of network:
%     N = size(G0,1);
%     spins = 1:N;
% 
%     heff = h;
%     Jeff = G0;
% 
%     Initialize things:
%     G = G0 ~= 0;
%     D = zeros(N-1,3);
% 
%     Find nodes with degree <= 2:
%     degs = sum(G);
%     inds = intersect(find(degs <= 2), find(degs >= 1));
%     inds(inds==keep(1)) = [];
%     inds(inds==keep(2)) = [];
%     inds(inds==keep(3)) = [];
% 
%     Loop until there are no more feasible nodes:
%     counter = 1;
% 
%     while ~isempty(inds) ~= 0
% 
%         Dtemp = zeros(1,3);
% 
%         Choose node to remove:
%         i = inds(randi(length(inds))); % Remove random node
%         i = inds(1); % Remove first node in list
%         i = inds(end); % Remove last node in list
%         Dtemp(1) = i;
% 
%         Neighbors of i:
%         Ni = spins(G(i,:));
%         Dtemp(2) = Ni(1);
% 
%         Remove node from network:
%         G(i,Ni) = 0;
%         G(Ni,i) = 0;
% 
%         heff(Ni(1)) = heff(Ni(1)) - log(exp(heff(i)) + 1) + log(exp(heff(i) + Jeff(i, Ni(1))) + 1);
% 
%         prevjkcon = 1;
%         If degree of i is 2 then connect two neighbors:
%         if length(Ni) == 2
%             Dtemp(3) = Ni(2);
% 
%             G(Ni(1), Ni(2)) = 1;
%             G(Ni(2), Ni(1)) = 1;
%             heff(Ni(2)) = heff(Ni(2)) - log(exp(heff(i)) + 1) + log(exp(heff(i) + Jeff(i, Ni(2))) + 1);
%             Jeff(Ni(1), Ni(2)) = Jeff(Ni(1),Ni(2)) + log(exp(heff(i)) + 1 ) - log( exp(Jeff(i, Ni(1)) + heff(i))+ 1)- log( exp(Jeff(i, Ni(2)) + heff(i))+ 1) + log( exp(Jeff(i, Ni(1)) + Jeff(i, Ni(2))+ heff(i))+ 1);
%             Jeff(Ni(2), Ni(1)) = Jeff(Ni(1),Ni(2));
%         end
% 
%         Remove node from network:
%         Jeff(i,Ni) = 0;
%         Jeff(Ni,i) = 0;
% 
%         Compute new feasible nodes to remove:
% 
%         degs(i) = degs(i) - 2;
%         if prevjkcon == 1
%             degs(Ni) = degs(Ni) - 1;
%         end
% 
%         binIntersect = (degs == 1) + (degs == 2);
% 
%         inds = spins(logical(binIntersect));
% 
% 
%         inds(inds==keep(1)) = [];
%         inds(inds==keep(2)) = [];
%         inds(inds==keep(3)) = [];
% 
%         D(counter,:) = Dtemp;
%         counter = counter + 1;
% 
%     end
% end
% 
% function [H] = Entropy(mean)
%     H = -mean.*log(mean)-(1-mean).*log(1-mean);
% end
% 
% function [spikes, probs] = Estimate_NumSpikeDist(J, h, kT)
%      USING 0,1 FOR SPINS
%     NEED TO ADD THERMALIZATION
% 
%     numSpins = length(h);
% 
%     initilize glauber with spin set with three up spins
%     spin_set = zeros(numSpins, 1);
%     spin_set(randperm(numel(spin_set), 3)) = 1;
% 
%     numIters = 2^12 ;
% 
%     probs = zeros(1,numel(spin_set)+1);
%     for iter = 1 : numIters
%         iter
%         for termal =  1:numel(spin_set)
%             Pick a random spin
%             Index = randi(numel(spin_set));
% 
%             Calculate energy change if this spin is flipped
%             dE = ((2*spin_set(Index)-1)*J(Index, :)* spin_set - J(Index, Index) + h(Index)*(2*spin_set(Index)-1));
% 
%             Boltzmann probability of flipping
%             prob = exp(-dE / kT)/(1+exp(-dE / kT));
% 
%             Spin flip condition
%             if rand() <= prob
%                 spin_set(Index) = 1 - spin_set(Index);
%             end
%         end
%         probs(sum(spin_set)+1) = probs(sum(spin_set)+1) + 1;
%     end
%     probs = probs/sum(probs);
%     spikes = 0:numel(spin_set);
% end
% 
% function deltaS = MI2(mean,corr)
% 
%     a = -mean.*log(mean) - (1-mean).*log(1-mean) - (mean.*log(mean)).' - ((1-mean).*log(1-mean)).';
%     b = log(corr.^corr); % Hack to make 0*log 0 = 0
%     c = (mean - corr).*log(mean - corr) + ((mean - corr).*log(mean - corr)).';
%     d = (1 - mean - mean.' +corr).*log((1 - mean - mean.' +corr));
%     deltaS = real(a + b + c + d);
% 
%     for i = 1:length(mean)
%         deltaS(i,i) = 0; %-mean(i)*log(mean(i)) - (1-mean(i))*log(1-mean(i));
%     end
% 
%     deltaS(isnan(deltaS)) = 0;
% end
% 
% function [mean,corr] = findStats(J,h)
%     Find the statistics of the GSP network defined in the Ising
%     interaction matrix J and external field vector h. Only finds
%     some of the correlation between nodes. Use correlations_GSP_01
%     instead for full correlation matrix. 
% 
% 
%     [child, parent1, parent2, heff, Jeff] = Decimate(J, h);
% 
%     child = flip(child,2);
%     parent1 = flip(parent1,2);
%     parent2 = flip(parent2,2); 
% 
%     c = 0 ;
%     for i = child
%         c = c + log(exp(heff(i))+1);
%     end
% 
%     mean = 1./(1+exp(-heff));
% 
%     corr = zeros(size(J));
% 
%     dim = size(child, 2);
% 
%     xi = child(2);
%     xj = parent1(2);
%     mean(xi) = (1-mean(xj))./(1+exp(-heff(xi))) + mean(xj)./(1+exp(-Jeff(xi, xj)-heff(xi)));
%     corr(xi,xj) = mean(xj)./(1+exp(-Jeff(xi, xj)-heff(xi)));
%     corr(xj,xi) = corr(xi,xj);
% 
%     for i = 3:dim
% 
%         xi = child(i);
%         xj = parent1(i);
%         xk = parent2(i);
%         if xk == 0
%             mean(xi) = (1-mean(xj))./(1+exp(-heff(xi))) + mean(xj)./(1+exp(-Jeff(xi, xj)-heff(xi)));
%             corr(xi,xj) = mean(xj)./(1+exp(-Jeff(xi, xj)-heff(xi)));
%             corr(xj,xi) = corr(xi,xj);
%         else
%             mean(xi) = (1-mean(xj) - mean(xk) + corr(xj,xk))./(1+exp(-heff(xi))) + (mean(xj) - corr(xj,xk))./(1+exp(-Jeff(xi, xj)-heff(xi))) + (mean(xk) - corr(xj,xk))./(1+exp(-Jeff(xi, xk)-heff(xi))) + (corr(xj,xk))./(1+exp(-Jeff(xi, xj)-Jeff(xi, xk)-heff(xi)));
%             corr(xi,xj) = (mean(xj) - corr(xj,xk))./(1+exp(-Jeff(xi, xj)-heff(xi))) + (corr(xj,xk))./(1+exp(-Jeff(xi, xj)-Jeff(xi, xk)-heff(xi)));
%             corr(xj,xi) = corr(xi,xj);
%             corr(xi,xk) = (mean(xk) - corr(xj,xk))./(1+exp(-Jeff(xi, xk)-heff(xi))) + (corr(xj,xk))./(1+exp(-Jeff(xi, xj)-Jeff(xi, xk)-heff(xi)));
%             corr(xk,xi) = corr(xi,xk);
%         end
% 
%     end
% 
% 
% end
% 
% function [child, parent1, parent2, heff, Jeff] = Decimate(J, h)
%     Decimate the GSP network defined in the connection matrix J
%     Not used in correlations_GSP_01.
%     numSpins = length(h);
%     G = graph(J, 'upper');
%     nodeNames = cellstr(num2str((1:numSpins).')); % Convert indices to cell array of strings
%     G.Nodes.Name = nodeNames; % Assign node names
%     Jeff = J;
% 
%     parent1 = [];
%     parent2 = [];
%     child = [];
%     nodenum = height(G.Nodes);
% 
%     heff = h;
%     while nodenum > 1
%         [m, Gind] = min(G.degree);
%         Gneigh = neighbors(G,Gind);
%         parents = outedges(G,Gind);
%         temp = G.Nodes.Name(Gind);
%         ind = str2num(temp{1}); % set as labled node
%         child = [child, ind];
% 
%         temp = G.Nodes.Name(Gneigh(1));
%         neigh(1) = str2num(temp{1}); % set as labled neighbor
%         heff(neigh(1)) = heff(neigh(1)) - log(exp(heff(ind)) + 1) + log(exp(heff(ind) + G.Edges.Weight(parents(1))) + 1);
% 
%         parent1 = [parent1, neigh(1)];
%         if m ==2
%             temp = G.Nodes.Name(Gneigh(2));
%             neigh(2) = str2num(temp{1}); % set as labled neighbor
%             heff(neigh(2)) = heff(neigh(2)) - log(exp(heff(ind)) + 1) + log(exp(heff(ind) + G.Edges.Weight(parents(2))) + 1);
%             parent2 = [parent2, neigh(2)];
%             edgeind = findedge(G,Gneigh(1),Gneigh(2));
%             if edgeind == 0
%                 G = addedge(G,Gneigh(1),Gneigh(2), log(exp(heff(ind)) + 1 ) - log( exp(G.Edges.Weight(parents(1)) + heff(ind))+ 1)- log( exp(G.Edges.Weight(parents(2)) + heff(ind))+ 1) + log( exp(G.Edges.Weight(parents(1))+G.Edges.Weight(parents(2)) + heff(ind))+ 1) ); 
%                 Jeff(neigh(1), neigh(2)) = log(exp(heff(ind)) + 1 ) - log( exp(G.Edges.Weight(parents(1)) + heff(ind))+ 1)- log( exp(G.Edges.Weight(parents(2)) + heff(ind))+ 1) + log( exp(G.Edges.Weight(parents(1))+G.Edges.Weight(parents(2)) + heff(ind))+ 1);
%                 Jeff(neigh(2), neigh(1)) = Jeff(neigh(1), neigh(2));
%             else
%                 weight = G.Edges.Weight(edgeind);
%                 Jeff(neigh(1), neigh(2)) = weight + log(exp(heff(ind)) + 1 ) - log( exp(G.Edges.Weight(parents(1)) + heff(ind))+ 1)- log( exp(G.Edges.Weight(parents(2)) + heff(ind))+ 1) + log( exp(G.Edges.Weight(parents(1))+G.Edges.Weight(parents(2)) + heff(ind))+ 1);
%                 Jeff(neigh(2), neigh(1)) = Jeff(neigh(1), neigh(2));
%                 G.Edges.Weight(edgeind) = weight + log(exp(heff(ind)) + 1 ) - log( exp(G.Edges.Weight(parents(1)) + heff(ind))+ 1)- log( exp(G.Edges.Weight(parents(2)) + heff(ind))+ 1) + log( exp(G.Edges.Weight(parents(1))+G.Edges.Weight(parents(2)) + heff(ind))+ 1) ;     
%             end
%         else
%             parent2 = [parent2, 0];
%         end
% 
%         remove node
%         G = rmnode(G, Gind);
%         nodenum = nodenum - 1;
% 
%     end
% 
%     [m, Gind] = min(G.degree);
%     temp = G.Nodes.Name(Gind);
%     ind = str2num(temp{1}); % set as labled node
%     child = [child, ind];
%     parent1 = [parent1,0];
%     parent2 = [parent2,0];
% 
% end
% 
% function [G, D] = decimate_GSP(G0)
%     Input: nxn unweighted adjacency matrix G0 representing a GSP network.
% 
%     Output: Randomly remove nodes of degree 1 or 2. If we remove a node of
%     degree 2 then we place a new edge between the two neighbors of the
%     removed node. Continue this process until no more nodes of degree 1 or 2
%     exist and return the final adjacency matrix G. Also return the order of
%     decimations in the numDec x 3 matrix D, where D(t,1) is the node removed
%     at step t and D(t,2) and D(t,3) are the two nodes from which the node was
%     removed. We note that D(t,3) = 0 if the removed node had degree 1.
% 
%     Size of network:
%     N = size(G0,1);
%     spins = 1:N;
% 
%     Initialize things:
%     G = G0 ~= 0;
%     D = zeros(N-1,3);
% 
%     Find nodes with degree <= 2:
%     degs = sum(G);
%     inds = intersect(find(degs <= 2), find(degs >= 1));
% 
%     Loop until there are no more feasible nodes:
%     counter = 1;
% 
%     while ~isempty(inds) ~= 0
% 
%         Dtemp = zeros(1,3);
% 
%         Choose node to remove:
%         randindx = randi(length(inds));
%         i = inds(randindx); % Remove random node
%         i = inds(1); % Remove first node in list
%         i = inds(end); % Remove last node in list
%         Dtemp(1) = i;
% 
% 
%         Neighbors of i:
%         Ni = spins(G(i,:));
%         Dtemp(2) = Ni(1);
% 
%         Remove node from network:
%         G(i,Ni) = 0;
%         G(Ni,i) = 0;
% 
% 
%         prevjkcon = 1;
%         If degree of i is 2 then connect two neighbors:
%         if length(Ni) == 2
%             Dtemp(3) = Ni(2);
% 
% 
%             G(Ni(1), Ni(2)) = 1;
%             G(Ni(2), Ni(1)) = 1;
%             prevjkcon = G(Ni(1), Ni(2));
% 
%         end
% 
%         Compute new feasible nodes to remove:
% 
%         degs(i) = degs(i) - 2;
%         if prevjkcon == 1
%             degs(Ni) = degs(Ni) - 1;
%         end
% 
%         binIntersect = (degs == 1) + (degs == 2);
% 
%         inds = spins(logical(binIntersect));
% 
%         D(counter,:) = Dtemp;
%         counter = counter + 1;
%         counter
% 
%     end
% end
% 
% function [Js, connum] =  Jpairs(J, D)
% 
% 
%     Jcon = J ~= 0;
%     Jcon = triu(Jcon);
% 
% 
%     Js = zeros(size(Jcon)); % I have no idea why i cant update Jcon
%     count = 0;
%     for iter = 1:length(D)
% 
%         i = D(iter,1);
%         j = D(iter,2);
%         k = D(iter,3);
% 
%         count = count + 1;
%         Js(i, j) = count;
%         if k ~=0 
%             count = count + 1;
%             Js(i, k) = count; 
%         end
% 
%     end
% 
%     Js = Js + Js.';
%     connum = count;
% 
% end
% 
% function [m, C, X, Z] = correlations_GSP_01(J, h)
%     Updated to be more memory efficient for all derivitives 11/3/23
%     Inputs: nxn matirx J and nx1 col vector h of external fields for a system
%     with 0,1 variables. We require that the interaction matrix J represents a
%     GSP network. We also take as input the order of node decimations (as
%     output by 'decimate_GSP'), where decimations(t,1) is the t^th node to be 
%     decimated and decimations(t,[2,3]) are its two neighbors. We note that
%     decimations(t,3) = 0 if the t^th node only has one neighbor.
% 
%     Output: nx1 vector m of magnetizations, calculated in two steps: (i)
%     decimate nodes of degree 1 or 2 down to a single node, and (ii)
%     recursivley calculate the magnetizations back out in reverse order. We
%     also compute the correlations C(i,j) between all nodes i and j and the
%     partition function Z. See Report 2 for details.
% 
%     NOTE: The difference between this function and 'correlations_GSP' is that
%     here consider systems with 0,1 variables.
% 
%     Number of nodes:
%     n = length(h);
% 
% 
%     disp(n)
%     Compute decimation order:
%     [Jdec, D] = decimate_GSP(J);
%     if sum(sum(Jdec)) ~= 0
%         error('Network was not able to be decimated!');
%     end
%     disp('decimated')
%     [Jcon, connum] =  Jpairs(J, D); % This function finds pairs that are directly connected in the network
% 
%     [child, parent1, parent2, h_eff, J_eff] = Decimate(J, h);
% 
% 
%     D = zeros(length(child)-1,3);
% 
%     D(:,1) = child(1:end-1);
%     D(:,2) = parent1(1:end-1);
%     D(:,3) = parent2(1:end-1);
% 
%     c = 0 ;
%     for i = child
%         c = c + log(exp(h_eff(i))+1);
%     end
% 
%     Initialize effective parameters:
%     J_eff = J;
%     h_eff = h;
%     c = 0;
%     disp('eff')
%     Loop through node decimations:
%     for t = 1:size(D,1)
% 
%         Node to decimate and its neighbors:
%         i = D(t,1);
%         j = D(t,2);
%         k = D(t,3);
% 
%         Update effective parameters on first neighbor:
%         c = c + log(exp(h_eff(i)) + 1);
%         h_eff(j) = h_eff(j) - log(exp(h_eff(i)) + 1) + log(exp(J_eff(i,j) + h_eff(i)) + 1);
% 
%         If node i has two neighbors:
%         if k ~= 0
% 
%             h_eff(k) = h_eff(k) - log(exp(h_eff(i)) + 1) + log(exp(J_eff(i,k) + h_eff(i)) + 1);
%             J_eff(j,k) = J_eff(j,k) + log(exp(h_eff(i)) + 1) - log(exp(J_eff(i,j) + h_eff(i)) + 1)...
%                 - log(exp(J_eff(i,k) + h_eff(i)) + 1) + log(exp(J_eff(i,j) + J_eff(i,k) + h_eff(i)) + 1);
%             J_eff(k,j) = J_eff(j,k);
% 
%         end
%     end
% 
%     Things we are going to need to compute:
%     m = zeros(n,1); % Magnetizations
%     C = zeros(n); % Correlations between nodes that interact
%     dm_dh = zeros(n); % Susceptibilities
%     dC_dh = zeros(connum,n); % Derivatives of correlations with respect to external fields
%     dh_dh = eye(n); % Derivatives of external fields with respect to external fields
%     dJ_dh = zeros(connum,n); % Derivatives of interactions with respect to external fields
%     dh_dJ = zeros(n,connum); % Derivatives of external fields with respect to interactions
%     Jpair = zeros(connum);
%     dJ_dJ = zeros(n,n,n,n); % Derivatives of interactions with respect to interactions
% 
%     Make self-derivatives unity:
%     for ind1 = 1:n
%         ind1
%         for ind2 = 1:n
%             if Jcon(ind1,ind2) ~= 0
%                 Jpair(Jcon(ind1,ind2),Jcon(ind1,ind2)) = 1;
%                 Jpair(Jcon(ind2,ind1),Jcon(ind1,ind2)) = 1;
%                 Jpair(Jcon(ind1,ind2),Jcon(ind2,ind1)) = 1;
%                 Jpair(Jcon(ind2,ind1),Jcon(ind2,ind1)) = 1;
%             end
%             dJ_dJ(ind1,ind2,ind1,ind2) = 1;
%             dJ_dJ(ind1,ind2,ind2,ind1) = 1;
%             dJ_dJ(ind2,ind1,ind1,ind2) = 1;
%             dJ_dJ(ind2,ind1,ind2,ind1) = 1;
%         end
%     end
% 
%     Compute things for final node:
%     i0 = j ;%child(end);
%     m(i0) = 1/(1 + exp(-h_eff(i0)));
%     C(i0,i0) = m(i0);
%     dm_dh(i0,i0) = exp(-h_eff(i0))/(1 + exp(-h_eff(i0)))^2;
% 
%     Compute partition function:
%     Z = exp(c)*(exp(h_eff(i0)) + 1);
% 
%     Loop over nodes in reverse order from which they were decimated:
%     for t_j = size(D,1):-1:1
%         t_j
%         Node to take derivative with respect to and its decimation neighbors:
%         j1 = D(t_j,1);
%         j2 = D(t_j,2);
%         j3 = D(t_j,3);
% 
%         Compute magnetizations and derivatives:
% 
%         If j1 only has one neighbor:
%         if j3 == 0
% 
%             Compute magnetizations and correlations:
%             m(j1) = (1 - m(j2))/(1 + exp(-h_eff(j1))) + m(j2)/(1 + exp(-J_eff(j1,j2) - h_eff(j1)));
%             C(j1,j1) = m(j1);
% 
%             C(j1,j2) = m(j2)/(1 + exp(-J_eff(j1,j2) - h_eff(j1)));
%             C(j2,j1) = C(j1,j2);
% 
%             Compute derivatives of h(j2) with respect to h(j1) and J(j1,j2):
%             dh_dh(j2,j1) = -1/(1 + exp(-h_eff(j1))) + 1/(1 + exp(-J_eff(j1,j2) - h_eff(j1)));
%             dh_dJ(j2,Jcon(j1,j2)) = 1/(1 + exp(-J_eff(j1,j2) - h_eff(j1)));
%             dh_dJ(j2,Jcon(j2,j1)) = dh_dJ(j2,Jcon(j1,j2));
% 
%             Compute derivative of external field of final node with respect
%             to h(j1) and J(j1,j2). Dependence goes through h(j2):
%             dh_dh(i0,j1) = dh_dh(i0,j2)*dh_dh(j2,j1);
%             dh_dJ(i0,Jcon(j1,j2)) = dh_dh(i0,j2)*dh_dJ(j2,Jcon(j1,j2));
%             dh_dJ(i0,Jcon(j2,j1)) = dh_dJ(i0,Jcon(j1,j2));
% 
%             Loop over nodes between final node and j1:
%             for t_i = size(D,1):-1:(t_j + 1)
% 
%                 Node to take derivative of and its decimation neighbors:
%                 i1 = D(t_i,1);
%                 i2 = D(t_i,2);
%                 i3 = D(t_i,3);
% 
%                 Compute derivatives of h(i1), J(i1,i2), and J(i1,i3) with
%                 respect to h(j1) and J(j1,j2). All dependencies go through h(j2):
%                 dh_dh(i1,j1) = dh_dh(i1,j2)*dh_dh(j2,j1);
%                 dh_dJ(i1,Jcon(j1,j2)) = dh_dh(i1,j2)*dh_dJ(j2,Jcon(j1,j2));
%                 dh_dJ(i1,Jcon(j2,j1)) = dh_dJ(i1,Jcon(j1,j2));
% 
%                 dJ_dh(Jcon(i1,i2),j1) = dJ_dh(Jcon(i1,i2),j2)*dh_dh(j2,j1);
%                 dJ_dh(Jcon(i2,i1),j1) = dJ_dh(Jcon(i1,i2),j1);
% 
% 
%                 Jpair(Jcon(i1,i2),Jcon(j1,j2)) = dJ_dh(Jcon(i1,i2),j2)*dh_dJ(j2,Jcon(j1,j2));
%                 Jpair(Jcon(i2,i1),Jcon(j1,j2)) = Jpair(Jcon(i1,i2),Jcon(j1,j2));
%                 Jpair(Jcon(i1,i2),Jcon(j2,j1)) = Jpair(Jcon(i1,i2),Jcon(j1,j2));
%                 Jpair(Jcon(i2,i1),Jcon(j2,j1)) = Jpair(Jcon(i1,i2),Jcon(j1,j2));
% 
%                 dJ_dJ(i1,i2,j1,j2) = dJ_dh(i1,i2,j2)*dh_dJ(j2,j1,j2);
%                 dJ_dJ(i2,i1,j1,j2) = dJ_dJ(i1,i2,j1,j2);
%                 dJ_dJ(i1,i2,j2,j1) = dJ_dJ(i1,i2,j1,j2);
%                 dJ_dJ(i2,i1,j2,j1) = dJ_dJ(i1,i2,j1,j2);
% 
%                 If i1 has two neighbors:
%                 if i3 ~= 0
% 
%                     dJ_dh(Jcon(i1,i3),j1) = dJ_dh(Jcon(i1,i3),j2)*dh_dh(j2,j1);
%                     dJ_dh(Jcon(i3,i1),j1) = dJ_dh(Jcon(i1,i3),j1);
% 
%                     Jpair(Jcon(i1,i3),Jcon(j1,j2)) = dJ_dh(Jcon(i1,i3),j2)*dh_dJ(j2,Jcon(j1,j2));
%                     Jpair(Jcon(i3,i1),Jcon(j1,j2)) = Jpair(Jcon(i1,i3),Jcon(j1,j2));
%                     Jpair(Jcon(i1,i3),Jcon(j2,j1)) = Jpair(Jcon(i1,i3),Jcon(j1,j2));
%                     Jpair(Jcon(i3,i1),Jcon(j2,j1)) = Jpair(Jcon(i1,i3),Jcon(j1,j2));
% 
%                     dJ_dJ(i1,i3,j1,j2) = dJ_dh(i1,i3,j2)*dh_dJ(j2,j1,j2);
%                     dJ_dJ(i3,i1,j1,j2) = dJ_dJ(i1,i3,j1,j2);
%                     dJ_dJ(i1,i3,j2,j1) = dJ_dJ(i1,i3,j1,j2);
%                     dJ_dJ(i3,i1,j2,j1) = dJ_dJ(i1,i3,j1,j2);
% 
%                 end
% 
%             end
% 
%         If j1 has two neighbors:
%         else
% 
%             Compute magnetizations and correlations:
%             m(j1) = (1 - m(j2) - m(j3) + C(j2,j3))/(1 + exp(-h_eff(j1)))...
%                 + (m(j2) - C(j2,j3))/(1 + exp(-J_eff(j1,j2) - h_eff(j1)))...
%                 + (m(j3) - C(j2,j3))/(1 + exp(-J_eff(j1,j3) - h_eff(j1)))...
%                 + C(j2,j3)/(1 + exp(-J_eff(j1,j2) - J_eff(j1,j3) - h_eff(j1)));
%             C(j1,j1) = m(j1);
% 
%             C(j1,j2) = (m(j2) - C(j2,j3))/(1 + exp(-J_eff(j1,j2) - h_eff(j1)))...
%                 + C(j2,j3)/(1 + exp(-J_eff(j1,j2) - J_eff(j1,j3) - h_eff(j1)));
%             C(j2,j1) = C(j1,j2);
% 
%             C(j1,j3) = (m(j3) - C(j2,j3))/(1 + exp(-J_eff(j1,j3) - h_eff(j1)))...
%                 + C(j2,j3)/(1 + exp(-J_eff(j1,j2) - J_eff(j1,j3) - h_eff(j1)));
%             C(j3,j1) = C(j1,j3);
% 
%             Compute derivatives of h(j2), h(j3), and J(j2,j3) with respect to
%             h(j1), J(j1,j2), and J(j1,j3):
%             dh_dh(j2,j1) = -1/(1 + exp(-h_eff(j1))) + 1/(1 + exp(-J_eff(j1,j2) - h_eff(j1)));
%             dh_dJ(j2,Jcon(j1,j2)) = 1/(1 + exp(-J_eff(j1,j2) - h_eff(j1)));
%             dh_dJ(j2,Jcon(j2,j1)) = dh_dJ(j2,Jcon(j1,j2));
% 
%             dh_dh(j3,j1) = -1/(1 + exp(-h_eff(j1))) + 1/(1 + exp(-J_eff(j1,j3) - h_eff(j1)));
%             dh_dJ(j3,Jcon(j1,j3)) = 1/(1 + exp(-J_eff(j1,j3) - h_eff(j1)));
%             dh_dJ(j3,Jcon(j3,j1)) = dh_dJ(j3,Jcon(j1,j3));
% 
%             dJ_dh(Jcon(j2,j3),j1) = 1/(1 + exp(-h_eff(j1))) - 1/(1 + exp(-J_eff(j1,j2) - h_eff(j1)))...
%                 - 1/(1 + exp(-J_eff(j1,j3) - h_eff(j1))) + 1/(1 + exp(-J_eff(j1,j2) - J_eff(j1,j3) - h_eff(j1)));
%             dJ_dh(Jcon(j3,j2),j1) = dJ_dh(Jcon(j2,j3),j1);
% 
% 
%             Jpair(Jcon(j2,j3),Jcon(j1,j2)) = -1/(1 + exp(-J_eff(j1,j2) - h_eff(j1)))...
%                 + 1/(1 + exp(-J_eff(j1,j2) - J_eff(j1,j3) - h_eff(j1)));
%             Jpair(Jcon(j2,j3),Jcon(j2,j1)) = Jpair(Jcon(j2,j3),Jcon(j1,j2));
%             Jpair(Jcon(j3,j2),Jcon(j1,j2)) = Jpair(Jcon(j2,j3),Jcon(j1,j2));
%             Jpair(Jcon(j3,j2),Jcon(j2,j1)) = Jpair(Jcon(j2,j3),Jcon(j1,j2));
% 
% 
%             dJ_dJ(j2,j3,j1,j2) = -1/(1 + exp(-J_eff(j1,j2) - h_eff(j1)))...
%                 + 1/(1 + exp(-J_eff(j1,j2) - J_eff(j1,j3) - h_eff(j1)));
%             dJ_dJ(j2,j3,j2,j1) = dJ_dJ(j2,j3,j1,j2);
%             dJ_dJ(j3,j2,j1,j2) = dJ_dJ(j2,j3,j1,j2);
%             dJ_dJ(j3,j2,j2,j1) = dJ_dJ(j2,j3,j1,j2);
% 
% 
%             Jpair(Jcon(j2,j3),Jcon(j1,j3)) = -1/(1 + exp(-J_eff(j1,j3) - h_eff(j1)))...
%                 + 1/(1 + exp(-J_eff(j1,j2) - J_eff(j1,j3) - h_eff(j1)));
%             Jpair(Jcon(j2,j3),Jcon(j3,j1)) = Jpair(Jcon(j2,j3),Jcon(j1,j3));
%             Jpair(Jcon(j3,j2),Jcon(j1,j3)) = Jpair(Jcon(j2,j3),Jcon(j1,j3));
%             Jpair(Jcon(j3,j2),Jcon(j3,j1)) = Jpair(Jcon(j2,j3),Jcon(j1,j3));
% 
%             dJ_dJ(j2,j3,j1,j3) = -1/(1 + exp(-J_eff(j1,j3) - h_eff(j1)))...
%                 + 1/(1 + exp(-J_eff(j1,j2) - J_eff(j1,j3) - h_eff(j1)));
%             dJ_dJ(j2,j3,j3,j1) = dJ_dJ(j2,j3,j1,j3);
%             dJ_dJ(j3,j2,j1,j3) = dJ_dJ(j2,j3,j1,j3);
%             dJ_dJ(j3,j2,j3,j1) = dJ_dJ(j2,j3,j1,j3);
% 
%             Compute derivative of external field of final node with respect
%             to h(j1), J(j1,j2), and J(j1,j3). Dependencies go through h(j2),
%             h(j3), and J(j2,j3):
%             dh_dh(i0,j1) = dh_dh(i0,j2)*dh_dh(j2,j1) + dh_dh(i0,j3)*dh_dh(j3,j1)...
%                 + dh_dJ(i0,Jcon(j2,j3))*dJ_dh(Jcon(j2,j3),j1);
%             dh_dJ(i0,Jcon(j1,j2)) = dh_dh(i0,j2)*dh_dJ(j2,Jcon(j1,j2)) + dh_dJ(i0,Jcon(j2,j3))*Jpair(Jcon(j2,j3),Jcon(j1,j2));  %dJ_dJ(j2,j3,j1,j2);
%             dh_dJ(i0,Jcon(j2,j1)) = dh_dJ(i0,Jcon(j1,j2));
%             dh_dJ(i0,Jcon(j1,j3)) = dh_dh(i0,j3)*dh_dJ(j3,Jcon(j1,j3)) + dh_dJ(i0,Jcon(j2,j3))*Jpair(Jcon(j2,j3),Jcon(j1,j3));  %dJ_dJ(j2,j3,j1,j3);
%             dh_dJ(i0,Jcon(j3,j1)) = dh_dJ(i0,Jcon(j1,j3));
% 
%             Loop over nodes between final node and j1:
%             for t_i = size(D,1):-1:(t_j + 1)
% 
%                 Node to take derivative of and its decimation neighbors:
%                 i1 = D(t_i,1);
%                 i2 = D(t_i,2);
%                 i3 = D(t_i,3);
% 
%                 Compute derivatives of h(i1), J(i1,i2), and J(i1,i3) with
%                 respect to h(j1), J(j1,j2) and J(j1,j3). All dependencies go
%                 through h(j2), h(j3), and J(j2,j3):
%                 dh_dh(i1,j1) = dh_dh(i1,j2)*dh_dh(j2,j1) + dh_dh(i1,j3)*dh_dh(j3,j1)...
%                     + dh_dJ(i1,Jcon(j2,j3))*dJ_dh(Jcon(j2,j3),j1);
%                 dh_dJ(i1,Jcon(j1,j2)) = dh_dh(i1,j2)*dh_dJ(j2,Jcon(j1,j2)) + dh_dJ(i1,Jcon(j2,j3))*Jpair(Jcon(j2,j3),Jcon(j1,j2)); %dJ_dJ(j2,j3,j1,j2);
%                 dh_dJ(i1,Jcon(j2,j1)) = dh_dJ(i1,Jcon(j1,j2));
%                 dh_dJ(i1,Jcon(j1,j3)) = dh_dh(i1,j3)*dh_dJ(j3,Jcon(j1,j3)) + dh_dJ(i1,Jcon(j2,j3))*Jpair(Jcon(j2,j3),Jcon(j1,j3)); %dJ_dJ(j2,j3,j1,j3);
%                 dh_dJ(i1,Jcon(j3,j1)) = dh_dJ(i1,Jcon(j1,j3));
% 
%                 dJ_dh(Jcon(i1,i2),j1) = dJ_dh(Jcon(i1,i2),j2)*dh_dh(j2,j1) + dJ_dh(Jcon(i1,i2),j3)*dh_dh(j3,j1)...
%                     + dJ_dh(Jcon(j2,j3),j1)*Jpair(Jcon(i1,i2),Jcon(j2,j3)); %dJ_dJ(i1,i2,j2,j3)*
%                 dJ_dh(Jcon(i2,i1),j1) = dJ_dh(Jcon(i1,i2),j1);
% 
%                 Jpair(Jcon(i1,i2),Jcon(j1,j2)) = dJ_dh(Jcon(i1,i2),j2)*dh_dJ(j2,Jcon(j1,j2)) + Jpair(Jcon(i1,i2),Jcon(j2,j3))*Jpair(Jcon(j2,j3),Jcon(j1,j2));
%                 Jpair(Jcon(i1,i2),Jcon(j2,j1)) = Jpair(Jcon(i1,i2),Jcon(j1,j2));
%                 Jpair(Jcon(i2,i1),Jcon(j1,j2)) = Jpair(Jcon(i1,i2),Jcon(j1,j2));
%                 Jpair(Jcon(i2,i1),Jcon(j2,j1)) = Jpair(Jcon(i1,i2),Jcon(j1,j2));
% 
%                 dJ_dJ(i1,i2,j1,j2) = dJ_dh(i1,i2,j2)*dh_dJ(j2,j1,j2) + dJ_dJ(i1,i2,j2,j3)*dJ_dJ(j2,j3,j1,j2);
%                 dJ_dJ(i1,i2,j2,j1) = dJ_dJ(i1,i2,j1,j2);
%                 dJ_dJ(i2,i1,j1,j2) = dJ_dJ(i1,i2,j1,j2);
%                 dJ_dJ(i2,i1,j2,j1) = dJ_dJ(i1,i2,j1,j2);
% 
%                 Jpair(Jcon(i1,i2),Jcon(j1,j3)) = dJ_dh(Jcon(i1,i2),j3)*dh_dJ(j3,Jcon(j1,j3)) + Jpair(Jcon(i1,i2),Jcon(j2,j3))*Jpair(Jcon(j2,j3),Jcon(j1,j3));
%                 Jpair(Jcon(i1,i2),Jcon(j3,j1)) = Jpair(Jcon(i1,i2),Jcon(j1,j3));
%                 Jpair(Jcon(i2,i1),Jcon(j1,j3)) = Jpair(Jcon(i1,i2),Jcon(j1,j3));
%                 Jpair(Jcon(i2,i1),Jcon(j3,j1)) = Jpair(Jcon(i1,i2),Jcon(j1,j3));
% 
%                 dJ_dJ(i1,i2,j1,j3) = dJ_dh(i1,i2,j3)*dh_dJ(j3,j1,j3) + dJ_dJ(i1,i2,j2,j3)*dJ_dJ(j2,j3,j1,j3);
%                 dJ_dJ(i1,i2,j3,j1) = dJ_dJ(i1,i2,j1,j3);
%                 dJ_dJ(i2,i1,j1,j3) = dJ_dJ(i1,i2,j1,j3);
%                 dJ_dJ(i2,i1,j3,j1) = dJ_dJ(i1,i2,j1,j3);
% 
%                 If i1 has two neighbors:
%                 if i3 ~= 0
% 
%                     dJ_dh(Jcon(i1,i3),j1) = dJ_dh(Jcon(i1,i3),j2)*dh_dh(j2,j1) + dJ_dh(Jcon(i1,i3),j3)*dh_dh(j3,j1)...
%                         + dJ_dh(Jcon(j2,j3),j1)*Jpair(Jcon(i1,i3),Jcon(j2,j3)); %dJ_dJ(i1,i3,j2,j3)*
%                     dJ_dh(Jcon(i3,i1),j1) = dJ_dh(Jcon(i1,i3),j1);
% 
%                     Jpair(Jcon(i1,i3),Jcon(j1,j2)) = dJ_dh(Jcon(i1,i3),j2)*dh_dJ(j2,Jcon(j1,j2)) + Jpair(Jcon(i1,i3),Jcon(j2,j3))*Jpair(Jcon(j2,j3),Jcon(j1,j2));
%                     Jpair(Jcon(i1,i3),Jcon(j2,j1)) = Jpair(Jcon(i1,i3),Jcon(j1,j2));
%                     Jpair(Jcon(i3,i1),Jcon(j1,j2)) = Jpair(Jcon(i1,i3),Jcon(j1,j2));
%                     Jpair(Jcon(i3,i1),Jcon(j2,j1)) = Jpair(Jcon(i1,i3),Jcon(j1,j2));
% 
%                     dJ_dJ(i1,i3,j1,j2) = dJ_dh(i1,i3,j2)*dh_dJ(j2,j1,j2) + dJ_dJ(i1,i3,j2,j3)*dJ_dJ(j2,j3,j1,j2);
%                     dJ_dJ(i1,i3,j2,j1) = dJ_dJ(i1,i3,j1,j2);
%                     dJ_dJ(i3,i1,j1,j2) = dJ_dJ(i1,i3,j1,j2);
%                     dJ_dJ(i3,i1,j2,j1) = dJ_dJ(i1,i3,j1,j2);
% 
% 
%                     Jpair(Jcon(i1,i3),Jcon(j1,j3)) = dJ_dh(Jcon(i1,i3),j3)*dh_dJ(j3,Jcon(j1,j3)) + Jpair(Jcon(i1,i3),Jcon(j2,j3))*Jpair(Jcon(j2,j3),Jcon(j1,j3));
%                     Jpair(Jcon(i1,i3),Jcon(j3,j1)) = Jpair(Jcon(i1,i3),Jcon(j1,j3));
%                     Jpair(Jcon(i3,i1),Jcon(j1,j3)) = Jpair(Jcon(i1,i3),Jcon(j1,j3));
%                     Jpair(Jcon(i3,i1),Jcon(j3,j1)) = Jpair(Jcon(i1,i3),Jcon(j1,j3));
% 
%                     dJ_dJ(i1,i3,j1,j3) = dJ_dh(i1,i3,j3)*dh_dJ(j3,j1,j3) + dJ_dJ(i1,i3,j2,j3)*dJ_dJ(j2,j3,j1,j3);
%                     dJ_dJ(i1,i3,j3,j1) = dJ_dJ(i1,i3,j1,j3);
%                     dJ_dJ(i3,i1,j1,j3) = dJ_dJ(i1,i3,j1,j3);
%                     dJ_dJ(i3,i1,j3,j1) = dJ_dJ(i1,i3,j1,j3);
% 
%                 end
% 
%             end
% 
%         end
% 
%         Compute susceptibilities:
% 
%         Compute susceptibility with final node (derivative of m(i0) with
%         respect to h(j1)). Dependence goes through h(i0):
%         dm_dh(i0,j1) = dm_dh(i0,i0)*dh_dh(i0,j1);
%         dm_dh(j1,i0) = dm_dh(i0,j1);
% 
%         Loop over nodes between final node and j1:
%         for t_i = size(D,1):-1:(t_j + 1)
% 
%             Node to take derivative of and its decimation neighbors:
%             i1 = D(t_i,1);
%             i2 = D(t_i,2);
%             i3 = D(t_i,3);
% 
%             If i1 only has one neighbor:
%             if i3 == 0
% 
%                 Useful quantities:
%                 l_h = 1/(1 + exp(-h_eff(i1)));
%                 l_hJ = 1/(1 + exp(-J_eff(i1,i2) - h_eff(i1)));
%                 dl_h = exp(-h_eff(i1))/(1 + exp(-h_eff(i1)))^2;
%                 dl_hJ = exp(-J_eff(i1,i2) - h_eff(i1))/(1 + exp(-J_eff(i1,i2) - h_eff(i1)))^2;
% 
%                 Compute susceptibility (derivative of m(i1) with respect to h(j1)):
%                 dm_dh(i1,j1) = dl_h*(1 - m(i2))*dh_dh(i1,j1) + dl_hJ*m(i2)*(dh_dh(i1,j1) + dJ_dh(Jcon(i1,i2),j1))...
%                     + (-l_h + l_hJ)*dm_dh(i2,j1);
%                 dm_dh(j1,i1) = dm_dh(i1,j1);
% 
%                 Compute derivative of C(i1,i2) with respect to h(j1):
%                 dC_dh(Jcon(i1,i2),j1) = dl_hJ*m(i2)*(dh_dh(i1,j1) + dJ_dh(Jcon(i1,i2),j1)) + l_hJ*dm_dh(i2,j1);
%                 dC_dh(Jcon(i2,i1),j1) = dC_dh(Jcon(i1,i2),j1);
% 
%             If i1 has two neighbors:
%             else
% 
%                 Useful quantities:
%                 l_h = 1/(1 + exp(-h_eff(i1)));
%                 l_hJ2 = 1/(1 + exp(-J_eff(i1,i2) - h_eff(i1)));
%                 l_hJ3 = 1/(1 + exp(-J_eff(i1,i3) - h_eff(i1)));
%                 l_hJJ = 1/(1 + exp(-J_eff(i1,i2) - J_eff(i1,i3) - h_eff(i1)));
%                 dl_h = exp(-h_eff(i1))/(1 + exp(-h_eff(i1)))^2;
%                 dl_hJ2 = exp(-J_eff(i1,i2) - h_eff(i1))/(1 + exp(-J_eff(i1,i2) - h_eff(i1)))^2;
%                 dl_hJ3 = exp(-J_eff(i1,i3) - h_eff(i1))/(1 + exp(-J_eff(i1,i3) - h_eff(i1)))^2;
%                 dl_hJJ = exp(-J_eff(i1,i2) - J_eff(i1,i3) - h_eff(i1))/(1 + exp(-J_eff(i1,i2) - J_eff(i1,i3) - h_eff(i1)))^2;
% 
%                 Compute susceptibility (derivative of m(i1) with respect to h(j1)):
%                 dm_dh(i1,j1) = dl_h*(1 - m(i2) - m(i3) + C(i2,i3))*dh_dh(i1,j1)...
%                     + dl_hJ2*(m(i2) - C(i2,i3))*(dh_dh(i1,j1) + dJ_dh(Jcon(i1,i2),j1))...
%                     + dl_hJ3*(m(i3) - C(i2,i3))*(dh_dh(i1,j1) + dJ_dh(Jcon(i1,i3),j1))...
%                     + dl_hJJ*C(i2,i3)*(dh_dh(i1,j1) + dJ_dh(Jcon(i1,i2),j1) + dJ_dh(Jcon(i1,i3),j1))...
%                     + (-l_h + l_hJ2)*dm_dh(i2,j1) + (-l_h + l_hJ3)*dm_dh(i3,j1)...
%                     + (l_h - l_hJ2 - l_hJ3 + l_hJJ)*dC_dh(Jcon(i2,i3),j1);
%                 dm_dh(j1,i1) = dm_dh(i1,j1);
% 
%                 Compute derivative of C(i1,i2) with respect to h(j1):
%                 dC_dh(Jcon(i1,i2),j1) = dl_hJ2*(m(i2) - C(i2,i3))*(dh_dh(i1,j1) + dJ_dh(Jcon(i1,i2),j1))...
%                     + dl_hJJ*C(i2,i3)*(dh_dh(i1,j1) + dJ_dh(Jcon(i1,i2),j1) + dJ_dh(Jcon(i1,i3),j1))...
%                     + l_hJ2*dm_dh(i2,j1) + (-l_hJ2 + l_hJJ)*dC_dh(Jcon(i2,i3),j1);
%                 dC_dh(Jcon(i2,i1),j1) = dC_dh(Jcon(i1,i2),j1);
% 
%                 Compute derivative of C(i1,i3) with respect to h(j1):
%                 dC_dh(Jcon(i1,i3),j1) = dl_hJ3*(m(i3) - C(i2,i3))*(dh_dh(i1,j1) + dJ_dh(Jcon(i1,i3),j1))...
%                     + dl_hJJ*C(i2,i3)*(dh_dh(i1,j1) + dJ_dh(Jcon(i1,i2),j1) + dJ_dh(Jcon(i1,i3),j1))...
%                     + l_hJ3*dm_dh(i3,j1) + (-l_hJ3 + l_hJJ)*dC_dh(Jcon(i2,i3),j1);
%                 dC_dh(Jcon(i3,i1),j1) = dC_dh(Jcon(i1,i3),j1);
% 
%             end
% 
%         end
% 
%     end
% 
%     Compute correlations between all nodes:
%     C = dm_dh + m*m';
%     C(logical(eye(n))) = m;
% 
%     Susceptibility:
%     X = C - m*m';
% end
% 
% function [meanspin, corr,threepoint] = Exact_Ising(J, h, kT)
%      USING 0,1 FOR SPINS
% 
%     numSpins = length(h);
% 
% 
%     z = 0; % Partition Function
%     probs = zeros(2^numSpins, 1);
%     for i = 0:2^numSpins -1
% 
%         spin_set = decimalToBinaryVector(i, numSpins).';
%         z = z + exp(-energyIsing(spin_set, J, h)/kT);
% 
%     end 
% 
%     meanspin = zeros(numSpins, 1);
%     corr = zeros(numSpins, numSpins);
%     threepoint = zeros(numSpins, numSpins,numSpins);
%     for i = 1:2^numSpins
%         spin_set = decimalToBinaryVector(i -1, numSpins).';
% 
%         probs(i) =  exp(-energyIsing(spin_set, J, h)/kT)/z;
%         meanspin = meanspin + spin_set.*probs(i);
%         corr =  corr + probs(i)*(spin_set)*spin_set.';
% 
% 
%         Create a 3D grid of matrices from the input vectors
%         [AA, BB, CC] = meshgrid(spin_set, spin_set, spin_set);
% 
%         Multiply the corresponding elements element-wise
%         tensor = AA .* BB .* CC;
% 
% 
%         threepoint = threepoint + probs(i)*tensor;
%     end 
% 
% 
% end
% 
% function Emean = energyIsing(spin, J, h)
%     %ENERGYISING Mean energy per spin.
%     %   Emean = ENERGYISING(spin, J) returns the mean energy per spin of the
%     %   configuration |spin|. |spin| is a matrix of +/- 1's. |J| is a scalar.
% 
%     Emean = -(1./2).*(spin.' * J * spin) + trace(J)  - h.'*spin;
% end
