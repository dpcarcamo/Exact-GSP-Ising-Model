%%
clear
clc

currentPath = pwd;
folders = split(currentPath, '\');
newPath = join(folders(1:length(folders)-2),"\");

addpath(strcat(newPath{1} , '\Model Data'))




% Add Helper Function to Path
newPath = join(folders(1:length(folders)-1),"\");
addpath(strcat(newPath{1}, '\Helper Function'))
load("natimg2800_M170717_MP033_2017-08-20_fitted.mat")
%%
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
histogram(Jtemp(Jtemp ~= 0),30,'Normalization','pdf', EdgeColor='none');
hold on
Jtemp = JGSP;
Jtemp(simul == 0) = 0;
histogram(Jtemp(Jtemp ~= 0),30, 'Normalization','pdf', EdgeColor='none');
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
[Y, E] = discretize(MIS2(MIS2>0),numbins);


[C,ia,ic] = unique(Y);
a_counts = accumarray(ic,1);

figure
loglog(E(Y(ia))/log(2), a_counts/sum(a_counts), '.', Color='k', MarkerSize=14)
xlabel('Mutual Information I_{ij} (bits)')
ylabel('Probability')
set(gca,'box','off') 

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

%% Accuracy of Model

[modelmean, modelcorr] = correlations_GSP_01(JGSP,hGSP);

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

bJ = JGSP ~= 0;
bJ = double(bJ);
metric = distances(graph(bJ));
metric(metric == 0) = nan;
metric(metric == Inf) = nan;

figure
plot(corr_coef(JGSP ~= 0 ),model_coef(JGSP ~= 0), '.')
hold on
plot(corr_coef(metric ==  2),model_coef(metric == 2), '.')
plot(linspace(-0.1,0.9),linspace(-0.1,0.9),'--', Color='black')
hold off
xlabel('Data Corr Coef', 'FontSize',16)
ylabel('Model Corr Coef', 'FontSize',16)
legend('Top. dis. 1','Top. dis. 2','Location','northwest', 'FontSize',14)

%% Plots for 1,2,3 topological distance


plotvariance(corr_coef(JGSP ~= 0 ),model_coef(JGSP ~= 0), 100, 'blue')
hold on
plotvariance(corr_coef(triu(metric==  2)), model_coef(triu(metric == 2)), 11 , 'red')
plotvariance(corr_coef(triu(metric ==  3)), model_coef(triu(metric == 3)), 6, [0.9290 0.6940 0.1250])
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

bJ = JTree ~= 0;
bJ = double(bJ);
metric = distances(graph(bJ));
metric(metric == 0) = nan;
metric(metric == Inf) = nan;

figure
plot(corr_coef(JTree ~= 0 ),model_coef(JTree ~= 0), '.')
hold on
plot(corr_coef(metric ==  2),model_coef(metric == 2), '.')
plot(corr_coef(metric ==  3),model_coef(metric == 3), '.')
plot(linspace(-0.1,0.9),linspace(-0.1,0.9),'--', Color='black')
hold off
xlabel('Data Corr Coef', 'FontSize',16)
ylabel('Model Corr Coef', 'FontSize',16)
legend('Top. dis. 1','Top. dis. 2','Top. dis. 3','Location','northwest', 'FontSize',14)

% Plots for 1,2,3 topological distance

figure
plotvariance(corr_coef(JTree ~= 0 ),model_coef(JTree ~= 0), 100, 'blue')
hold on
plotvariance(corr_coef(triu(metric ==  2)), model_coef(triu(metric == 2)), 11 , 'red')
plotvariance(corr_coef(triu(metric ==  3)), model_coef(triu(metric == 3)), 6, [0.9290 0.6940 0.1250])
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

bJ = JGSPdist ~= 0;
bJ = double(bJ);
metric = distances(graph(bJ));
metric(metric == 0) = nan;
metric(metric == Inf) = nan;

figure
plot(corr_coef(JGSPdist ~= 0 ),model_coef(JGSPdist ~= 0), '.')
hold on
plot(corr_coef(metric ==  2),model_coef(metric == 2), '.')
plot(corr_coef(metric ==  3),model_coef(metric == 3), '.')
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
plotvariance(corr_coef(triu(metric ==  2)), model_coef(triu(metric == 2)), 11 , 'red')
plotvariance(corr_coef(triu(metric ==  3)), model_coef(triu(metric== 3)), 6, [0.9290 0.6940 0.1250])
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


bJ = JTreedist ~= 0;
bJ = double(bJ);
metric = distances(graph(bJ));
metric(metric == 0) = nan;
metric(metric == Inf) = nan;
figure
plot(corr_coef(JTreedist ~= 0 ),model_coef(JTreedist ~= 0), '.')
hold on
plot(corr_coef(metric ==  2),model_coef(metric == 2), '.')
plot(corr_coef(metric ==  3),model_coef(metric == 3), '.')
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
plotvariance(corr_coef(triu(metric ==  2)), model_coef(triu(metric == 2)), 11 , 'red')
plotvariance(corr_coef(triu(metric ==  3)), model_coef(triu(metric == 3)), 6, [0.9290 0.6940 0.1250])
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

bJ = JGSPrand ~= 0;
bJ = double(bJ);
metric = distances(graph(bJ));
metric(metric == 0) = nan;
metric(metric == Inf) = nan;
figure

figure
plot(corr_coef(JGSPrand ~= 0 ),model_coef(JGSPrand ~= 0), '.')
hold on
plot(corr_coef(metric ==  2),model_coef(metric == 2), '.')
plot(corr_coef(metric ==  3),model_coef(metric == 3), '.')
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
plotvariance(corr_coef(triu(metric ==  2)), model_coef(triu(metric == 2)), 11 , 'red')
plotvariance(corr_coef(triu(metric ==  3)), model_coef(triu(metric == 3)), 6, [0.9290 0.6940 0.1250])
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

%% Triplets Accuracy

tic
[G,D] = decimate_GSP(JGSP);
toc

Data = zeros(2,100);
count = 1;
for keep = D.'
    
    tic
    [Gtest, Dtest, Jtest, htest] = decimate_GSP_carefully(JGSP,hGSP, keep);
    toc

    [meanspin, corr,threepoint] = Exact_Ising(Jtest(keep,keep), htest(keep).', 1);
    
    Data(1,count) = threepoint(1,2,3);
    
    Data(2,count) = (binary_M(keep(1),:).*binary_M(keep(2),:))*binary_M(keep(3),:).'/size(binary_M,2);

    Data(1,count)
    Data(2,count)
    count = count + 1

    
end

%%
plot(Data(:,1),Data(:,2), '.')
hold on
plot([min(Data(:,1)),max(Data(:,1))], [min(Data(:,1)),max(Data(:,1))], '--',Color='k')
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


%% Prob Active v. Effective Field


hGSPeffmatrix = hGSP + binary*JGSP;

[uniqueValues, uniqueIndices] = findUniqueEntries(hGSPeffmatrix);
%%
hGSPeff = zeros(20,1);
Pongivenheff = zeros(20,1);

count = 1;
for i = 1:num_nuerons  
    
    vals = uniqueValues{i};
    indexes = uniqueIndices{i};

    for j = 1:length(vals)

        hGSPeff(count) = vals(j);
        num_heff = sum(indexes == j);
        num_heff_on = sum(binary(indexes == j,i)); 

        Pongivenheff(count) = num_heff_on/num_heff;
        count = count + 1;

    end

end

%%

figure
[~, I] = sort(hGSPeff);
%plot(hGSPeff(I),datamean(I), '.')
plot(hGSPeff,Pongivenheff,'.')
hold on
plot(hGSPeff(I),1./(1+exp(-hGSPeff(I))))
hold off
xlabel('heff')
ylabel('Prob On')

%% Averaging 

[C,ia,ic] = unique(Pongivenheff);
avhGSPeff = zeros(3,1);

for i = 1:length(ia)
    avhGSPeff(i) = mean(hGSPeff(ic == i));
end
%plot(avhGSPeff,C,'.')
hold on
plot(hGSPeff(I),1./(1+exp(-hGSPeff(I))), 'LineWidth',1)
plotvariance(avhGSPeff,C,10,'red')
hold off
xlabel('heff')
ylabel('Prob On')
xlim([-8,8])
ylim([0,1])

%% Bining x axis

figure
[~, I] = sort(hGSPeff);
%plot(hGSPeff(I),datamean(I), '.')
plot(hGSPeff,Pongivenheff,'.')
hold on
%plot(hGSPeff(I),1./(1+exp(-hGSPeff(I))))
plotvariance(hGSPeff,Pongivenheff,15,'red')
hold off
xlabel('heff')
ylabel('Prob On')


%%
hTreeeff = hTree + datamean*JTree;

figure
[~, I] = sort(hTreeeff);
plot(hTreeeff(I),datamean(I), '.')
hold on
plot(hTreeeff(I),1./(1+exp(-hTreeeff(I))))
hold off

%% Entropy Drop Prediction

bJGSP = JGSP ~= 0;
bJGSP2 = bJGSP*bJGSP;
bJGSP2(bJGSP == 0) = 0;


Hind = Entropy(datamean);

mutual_info = MI2(datamean,datacorr_pseudo);

tic
[G,D] = decimate_GSP(JGSP);
toc

Data = zeros(2,100);
count = 1;

deltaH = 0;
for keep = D(1:size(D,1)-1,:).'
    
    % Both methods work
    %[Gtest, Dtest, Jtest, htest] = decimate_GSP_carefully(JGSP,hGSP, keep);
    [htest, Jtest, Htri, Htri2] = find_GSP_update(datamean(keep).',datacorr_pseudo(keep,keep));

    %Htri = Entropy_GSP(Jtest(keep,keep), htest(keep).', datamean(keep),datacorr_pseudo(keep,keep));

    deltaH = deltaH + sum(Hind(keep)) - Htri
    count = count + 1
    
end


twice_edge = mutual_info.*(bJGSP2-1);

twice_edge(bJGSP == 0) = 0;

deltaH = deltaH - sum(sum(twice_edge))/2
sum(Hind) - HGSP


function [uniqueValues, uniqueIndices] = findUniqueEntries(matrix)
    % Get the number of columns in the matrix
    numColumns = size(matrix, 2);

    % Initialize cell arrays to store unique values and their indices for each column
    uniqueValues = cell(1, numColumns);
    uniqueIndices = cell(1, numColumns);

    % Iterate through each column of the matrix
    for col = 1:numColumns
        % Find unique values and their indices in the current column
        [uniqueValues{col}, ~, uniqueIndices{col}] = unique(matrix(:, col));

    end
end

