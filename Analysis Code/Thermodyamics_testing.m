clear
clc
%%
currentPath = pwd;
folders = split(currentPath, '\');
newPath = join(folders(1:length(folders)-2),"\");

addpath(strcat(newPath{1} , '\Stringer Data\Data'))
%ExpData = matfile('natimg2800_M170717_MP033_2017-08-20.mat');

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


% Binarizing Data
binary = (Response - mu)> 2*sigma;
%imshow(binary)


% Filter Data
spiking_patterns = binary.';
num_bins = size(spiking_patterns,2);
num_nuerons = size(spiking_patterns,1);

datacorr = spiking_patterns*spiking_patterns.'/num_bins;
datamean = mean(spiking_patterns.');

datacorr_pseudo = datacorr + ones(num_nuerons)/(num_bins +1);


Hind = sum(Entropy(datamean));




%% Distance

distanceSMetric = squareform(pdist(Locations));
distanceMetric = distanceSMetric(5006,:);

%%


[sortDistances, I] = sort(distanceMetric);


%%

Ts = linspace(0.1,4);
Specifics = zeros(size(Ts,2),4,1000);
count = 1;
for j = unique([10,100,1000])
    i = count;
    for k = 1:round(num_nuerons/j)
        
        distanceMetric = distanceSMetric(randi(num_nuerons),:);
        [sortDistances, I] = sort(distanceMetric);
        region   = (distanceMetric <=  sortDistances(j));
        
        reduceddatamean = datamean(region);
        reduceddatacorr = datacorr_pseudo(region, region);
        reducedLocation = Locations(region,:);
        
        Hind = sum(Entropy(reduceddatamean));
            
        % Find Optimal GSP
        
        tic
        [hGSP, JGSP, HGSP, NewEnt2] = find_GSP_update(reduceddatamean, reduceddatacorr);
        toc

        %Entropy Drops
        
        
        
        count2 = 1;
        for T = Ts

            [Energy, Specific] = Thermodynamics(JGSP,hGSP, T);
            Specifics(count2, i,k) = Specific;
            count2 = count2 + 1;
        
        end
    end
    count = count + 1
    j
end
%%
nSpecifics = Specifics;
nSpecifics(Specifics == 0) = nan;

%%
nSpecifics = nanmean(nSpecifics, 3);

%%

currentPath = pwd;
folders = split(currentPath, '\');
newPath = join(folders(1:length(folders)-2),"\");

addpath(strcat(newPath{1} , '\Model Data'))
listing = dir(strcat(newPath{1} , '\Model Data'));

filename = listing(8).name

%
load(filename)

%%
Ts = linspace(0.1,4);
count2 = 1;
for T = Ts

    [Energy, Specific] = Thermodynamics(JGSP,hGSP, T);
    nSpecifics(count2, 4) = Specific;
    count2 = count2 + 1

end
%%
clc

Ns = [10,100,1000,10144];
alphas = [0.25,0.5,0.75,1];
cmap = colormap(winter);
for c = 1:4
    plot(Ts,nSpecifics(:,c)/Ns(c), 'Color',cmap(c*60,:),'LineWidth',2 )
    hold on
end
plot(ones(100,1),linspace(0,2.5), '--', 'Color', [0.6,0.6,0.6],'LineWidth',1)
hold off
xlabel('Temperature (k_bT)', 'FontSize',18)
ylabel('Specific Heat per Neuron (C(T)/N)', 'FontSize',18)
legend('10','100','1000','10144')
axis square


