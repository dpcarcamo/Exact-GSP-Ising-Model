clear
clc

addpath('Stringer Data\Data')
ExpData = matfile('natimg2800_M170717_MP033_2017-08-20.mat');


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


function [JGSPdist, hGSPdist, H] = Minimum_Distance_GSP(datamean,datacorr_pseudo,Locations)

    distanceMetric = squareform(pdist(Locations));
    
    numSpins = size(distanceMetric,1);
    spins = 1:numSpins;
    
    distanceMetric = distanceMetric + diag(inf*ones(numSpins,1));
    
    JGSPdist = zeros(numSpins,numSpins);
    
    hGSPdist = log(datamean./(1-datamean));
    
    dims = size(distanceMetric);
    
    deltaS = 1./distanceMetric;
    
    deltaS = deltaS - diag(diag(deltaS));
    
    [M,I] = max(deltaS, [], 'all', 'linear');
        
    [row, col] = ind2sub(dims,I);
    
    %Initialize first connection
    
    deltaS(row, col) = 0;
    deltaS(col, row) = 0;
    hGSPdist(row) = log( (datamean(row) -datacorr_pseudo(row,col))/ (1 + datacorr_pseudo(row,col) - datamean(row) -datamean(col)) );
    JGSPdist(row,col) = log( datacorr_pseudo(row,col)/ (datamean(col) - datacorr_pseudo(row,col))) - hGSPdist(row);
    JGSPdist(col, row) = JGSPdist(row, col);
    
    hGSPdist(col) = hGSPdist(col) + log(exp(hGSPdist(row)) + 1)- log(exp(JGSPdist(row,col) + hGSPdist(row)  )+ 1);
    connectedSpins = [row, col];
    
    count = 1
    
    Hind = sum(Entropy(datamean));
    H2 = MI2(datamean([row,col]),datacorr_pseudo([row,col],[row,col]));
    H = Hind - H2(1,2);         
    
    Addedpairs = zeros(2*numSpins - 3, 2); 
    Addedpairs(1,1) = col;
    Addedpairs(1,2) = row;
    distanceTri = zeros(numSpins, 2*numSpins - 3);
    
    spins(connectedSpins) = [];
    for spin = spins
        distanceTri(spin, 1) = 1./(distanceMetric(spin, row) + distanceMetric(spin, col));
    end
    
    
    for step = 1:numSpins -2
        
        [M,I] = max(distanceTri, [],"all", "linear");
    
        dim = size(distanceTri);
        [i,col] = ind2sub(dim,I);
    
        jk = Addedpairs(col,:);
        j = jk(1);
        k = jk(2);
        
    
        connectedSpins = [connectedSpins, i]; %Update list of already connected spins
    
        
        
        %-----%
    
        m = datamean([i,j,k]);
        C = datacorr_pseudo([i,j,k],[i,j,k]);
        
        step_size = 1;
        [h1, J12, J13] = inverse_Ising_GSP_01_helper(m, C, step_size);
    
        hGSPdist(i) = h1;
        Jij = J12;
        Jik = J13;
    
        %-----%
    
        JGSPdist(i, j) = Jij;
        JGSPdist(j, i) = Jij;
        JGSPdist(i, k) = Jik;
        JGSPdist(k, i) = Jik;
    
    
        hjold = hGSPdist(j);
        hkold = hGSPdist(k);
        Jjkold = JGSPdist(j,k);
    
        hj = hGSPdist(j) + log(exp(hGSPdist(i))+1) - log(exp(Jij + hGSPdist(i)) + 1);
        hk = hGSPdist(k) + log(exp(hGSPdist(i))+1) - log(exp(Jik + hGSPdist(i)) + 1);
        Jjk = JGSPdist(j,k) - log(exp(hGSPdist(i)) + 1 ) + log( exp(Jij + hGSPdist(i))+ 1) + log( exp(Jik + hGSPdist(i))+ 1) - log( exp(Jij + Jik + hGSPdist(i))+ 1) ;
    
    
        hGSPdist(j) = hj;
        hGSPdist(k) = hk;
        JGSPdist(j,k) = Jjk;
        JGSPdist(k,j) = JGSPdist(j,k);
    
        deltaH = -log(exp(hGSPdist(i))+1) + JGSPdist(i,j)*datacorr_pseudo(i,j) + JGSPdist(i,k)*datacorr_pseudo(i,k) + hGSPdist(i)*datamean(i) - datamean(i)*log(datamean(i))- (1-datamean(i))*log(1-datamean(i)) + (JGSPdist(j,k) - Jjkold)*datacorr_pseudo(j,k) + (hGSPdist(j) - hjold)*datamean(j) + (hGSPdist(k)- hkold)*datamean(k);
    
        pair = sort([i, j]);
        Addedpairs(2*step,1) = pair(1);
        Addedpairs(2*step,2) = pair(2);
        pair = sort([i, k]);
        Addedpairs(2*step+1,1) = pair(1);
        Addedpairs(2*step+1,2) = pair(2);
    
        spins = 1:numSpins;
        spins(connectedSpins) = [];
        for spin = spins
            distanceTri(spin, 2*step) = 1./(distanceMetric(spin, i) + distanceMetric(spin, j));
            distanceTri(spin, 2*step+1) = 1./(distanceMetric(spin, i) + distanceMetric(spin, k));
        end
    
        
        distanceTri(i, :) = 0;
        step
        H = H - deltaH  ;
    end
    H/Hind

end 

function [JTreedist, hTreedist, H] = Minimum_Distance_Tree(datamean,datacorr_pseudo,Locations)

    distanceMetric = squareform(pdist(Locations));
    
    numSpins = size(distanceMetric,1);
    spins = 1:numSpins;
    
    distanceMetric = distanceMetric - diag(inf*ones(numSpins,1));
    
    JTreedist = zeros(numSpins,numSpins);
    
    hTreedist = log(datamean./(1-datamean));
    
    dims = size(distanceMetric);
    
    deltaS = 1./distanceMetric;
    
    deltaS = deltaS - diag(diag(deltaS));
    
    [M,I] = max(deltaS, [], 'all', 'linear');
        
    [row, col] = ind2sub(dims,I);
    
    %Initialize first connection
    
    deltaS(row, col) = 0;
    deltaS(col, row) = 0;
    hTreedist(row) = log( (datamean(row) -datacorr_pseudo(row,col))/ (1 + datacorr_pseudo(row,col) - datamean(row) -datamean(col)) );
    JTreedist(row,col) = log( datacorr_pseudo(row,col)/ (datamean(col) - datacorr_pseudo(row,col))) - hTreedist(row);
    JTreedist(col, row) = JTreedist(row, col);
    
    hTreedist(col) = hTreedist(col) + log(exp(hTreedist(row)) + 1)- log(exp(JTreedist(row,col) + hTreedist(row)  )+ 1);
    numSpins = numSpins -2;
    connectedSpins = [row, col];
    spins = spins(spins~=row );
    spins = spins(spins~=col );
    
    count = 1
    
    Hind = sum(Entropy(datamean));
    H2 = MI2(datamean,datacorr_pseudo);
    H = Hind - H2(row,col);            
    
    for i= 1:numSpins
        [M, row, col, trow, tcol] = find_max_in_selected(deltaS, spins, connectedSpins); % Find maximum change in entropy
        connectedSpins = [connectedSpins, row]; %Update list of already connected spins
        spins(trow) = [];
    
        deltaS(row, col) = 0;
        deltaS(col, row) = 0;
        hTreedist(row) = log( (datamean(row) -datacorr_pseudo(row,col))/ (1 + datacorr_pseudo(row,col) - datamean(row) -datamean(col)) );
        JTreedist(row,col) = log( datacorr_pseudo(row,col)/ (datamean(col) - datacorr_pseudo(row,col))) - hTreedist(row);
        JTreedist(col, row) = JTreedist(row, col);
    
        hTreedist(col) = hTreedist(col) + log(exp(hTreedist(row)) + 1)- log(exp(JTreedist(row,col) + hTreedist(row)  )+ 1);
        count = count + 1
        H = H - H2(row,col) ; 
    end
    H/Hind

end

function [J, h, H] = RandomGSPFit(datamean, datacorr)
    % Function to fit a random GSP network onto data. 

    numSpins = length(datamean);


    J = zeros(numSpins); % Initialize random Spin matrix
    
    h = randn(numSpins,1);


    linidx = randi((numSpins-1)*(numSpins-2)/2);
    p = (sqrt(1+8*linidx)-1)/2;
    idx0 = floor(p);

    row = idx0 + 1 + 1;
    col = linidx - idx0*(idx0 + 1)/2 ;


    connectedSpins = [row, col];
    
    h(row) = log( (datamean(row) -datacorr(row,col))/ (1 + datacorr(row,col) - datamean(row) -datamean(col)) );
    J(row,col) = log( datacorr(row,col)/ (datamean(col) - datacorr(row,col))) - h(row);
    J(col, row) = J(row, col);
    
    h(col) = h(col) + log(exp(h(row)) + 1)- log(exp(J(row,col) + h(row)  )+ 1);


    spins = 1:numSpins;
    spins(connectedSpins) = [];

    Hind = sum(Entropy(datamean));
    H2 = MI2(datamean([row,col]),datacorr([row,col],[row,col]));
    H = Hind - H2(1,2);


    addedpairs = zeros(2*numSpins-3,2);
    addedpairs(1,1) = row;
    addedpairs(1,2) = col;

    for step = 3:numSpins

        randidx = randi(2*step-5);

        i = spins(randi(numSpins-step+1));

        j = addedpairs(randidx,1);
        k = addedpairs(randidx,2);

        addedpairs(2*step-4,1) = i;
        addedpairs(2*step-4,2) = j;
        addedpairs(2*step-3,1) = i;
        addedpairs(2*step-3,2) = k;
        
        connectedSpins = [connectedSpins, i]; %Update list of already connected spins
    
            
        %-----%
    
        m = datamean([i,j,k]);
        C = datacorr([i,j,k],[i,j,k]);
        
        step_size = 1;
        [h1, J12, J13] = inverse_Ising_GSP_01_helper(m, C, step_size);
    
        h(i) = h1;
        Jij = J12;
        Jik = J13;
    
        %-----%
    
        J(i, j) = Jij;
        J(j, i) = Jij;
        J(i, k) = Jik;
        J(k, i) = Jik;
    
    
        hjold = h(j);
        hkold = h(k);
        Jjkold = J(j,k);
    
        hj = h(j) + log(exp(h(i))+1) - log(exp(Jij + h(i)) + 1);
        hk = h(k) + log(exp(h(i))+1) - log(exp(Jik + h(i)) + 1);
        Jjk = J(j,k) - log(exp(h(i)) + 1 ) + log( exp(Jij + h(i))+ 1) + log( exp(Jik + h(i))+ 1) - log( exp(Jij + Jik + h(i))+ 1) ;
    
    
        h(j) = hj;
        h(k) = hk;
        J(j,k) = Jjk;
        J(k,j) = J(j,k);
    
        deltaH = -log(exp(h(i))+1) + J(i,j)*datacorr(i,j) + J(i,k)*datacorr(i,k) + h(i)*datamean(i) - datamean(i)*log(datamean(i))- (1-datamean(i))*log(1-datamean(i)) + (J(j,k) - Jjkold)*datacorr(j,k) + (h(j) - hjold)*datamean(j) + (h(k)- hkold)*datamean(k);
    

        spins = 1:numSpins;
        spins(connectedSpins) = [];
    

        step
        H = H - deltaH;  

    end
    

end

function [spikes, probs] = Estimate_NumSpikeDist(J, h, kT)
    %  USING 0,1 FOR SPINS
    % NEED TO ADD THERMALIZATION

    numSpins = length(h);

    % initilize glauber with spin set with three up spins
    spin_set = zeros(numSpins, 1);
    spin_set(randperm(numel(spin_set), 3)) = 1;

    numIters = 2^16 ;

    probs = zeros(1,numel(spin_set)+1);
    for iter = 1 : numIters
        for termal =  1:numel(spin_set)
            % Pick a random spin
            Index = randi(numel(spin_set));
            
            % Calculate energy change if this spin is flipped
            dE = ((2*spin_set(Index)-1)*J(Index, :)* spin_set - J(Index, Index) + h(Index)*(2*spin_set(Index)-1));
            
            % Boltzmann probability of flipping
            prob = exp(-dE / kT)/(1+exp(-dE / kT));
            
            % Spin flip condition
            if rand() <= prob
                spin_set(Index) = 1 - spin_set(Index);
            end
        end
        probs(sum(spin_set)+1) = probs(sum(spin_set)+1) + 1;
    end
    probs = probs/sum(probs);
    spikes = 0:numel(spin_set);
end

function H = Entropy_GSP(J,h, ExactMean, ExactCorr)
    % Currently Incorrect
  
    Jeff = J;
    heff = h;

    [G, D] = decimate_GSP(J);

    H = sum(Entropy(ExactMean)); %Entropy Independent model
    H2 = MI2(ExactMean,ExactCorr);

    
    c = 0 ;


    for step = 1:size(D,1)

        i = D(step, 1);
        j = D(step, 2);
        k = D(step, 3);
        
        heff(j) = heff(j) - log(exp(heff(i))+1) + log(exp(Jeff(j, i)+ heff(i))+1);

        if k ~= 0
            heff(k) = heff(k) - log(exp(heff(i))+1) + log(exp(Jeff(k, i)+ heff(i))+1);
            Jeff(j,k) = Jeff(j,k) + log(exp(heff(i))+1) - log(exp(Jeff(i,j)+heff(i))+1) - log(exp(Jeff(i,k)+heff(i))+1) + log(exp(Jeff(i,j)+ Jeff(i,k)+heff(i))+1);
            Jeff(k, j) = Jeff(j, k);
        end


    end

    D = flip(D);


    dim = size(D, 1);

    for i = 1:length(h)
        c = c + log(exp(heff(i))+1);
    end

    
    for i = 1:dim  
        xi = D(i,1);
        xj = D(i,2);
        xk = D(i,3);
        if xk == 0
            H = H - H2(xi,xj);
        else
            MI3 = -log(exp(heff(xi))+1) + Jeff(xi,xj)*ExactCorr(xi,xj)+Jeff(xi,xk)*ExactCorr(xi,xk) + heff(xi)*ExactMean(xi) - ExactMean(xi)*log(ExactMean(xi))-(1-ExactMean(xi))*log(1-ExactMean(xi))+(-log(exp(heff(xi))+1) + log(exp(Jeff(xi,xj)+heff(xi))+1) + log(exp(Jeff(xi,xk)+heff(xi))+1) - log(exp(Jeff(xi,xj)+ Jeff(xi,xk)+ heff(xi))+ 1))*ExactCorr(xj,xk) + (log(exp(heff(xi))+1) - log(exp(Jeff(xi,xj)+ heff(xi))+ 1))*ExactMean(xj)+ (log(exp(heff(xi))+1) - log(exp(Jeff(xi,xk)+ heff(xi))+ 1))*ExactMean(xk);
            H = H - MI3;
        end
    end

end

function [threepoint] = Chow_Liu(i, j, k, mean, corr)
    % Aproximate a triangle by its subtree that maximizes the mutual info

    spins = [i, j, k];
    MI = MI2(mean(spins),corr(spins,spins));
    MI = triu(MI);
    MI = MI + tril(ones(3));
    [M, I] = min(MI, [], 'all', 'linear');
    [row, col] = ind2sub(size(MI),I);
    MI(row, col) = 0;
    MI = triu(MI) + triu(MI).' - 2*diag(diag(ones(3)));

    [leaves, node] = Findleaves(MI); % Find the remaining tree

    [GC,GR] = groupcounts([leaves, node]');
    if length(GR) == 3

        for iter = 1:3
            if GC(iter) == 2
                first = spins(mod(GR(iter) -2,3) + 1);
                middle = spins(GR(iter));
                last = spins(mod(GR(iter), 3) + 1) ;
            end
        end
     
    
        prob1 = corr(first, middle)*corr(middle, last)/mean(middle)
    elseif length(GR) == 2
        first = spins(GC(1));
        last = spins(GC(2));
        if 1 ~= GC(1) || 1 ~= GC(2)
            prob1 = corr(first,last)*mean(spins(1));
        end
        if 2 ~= GC(1) || 2 ~= GC(2)
            prob1 = corr(first,last)*mean(spins(2));
        end
        if 3 ~= GC(1) || 3 ~= GC(2)
            prob1 = corr(first,last)*mean(spins(3));
        end
    else 
        prob1 = mean(spins(1))*mean(spins(2))*mean(spins(3));
    end
    threepoint = zeros(3,3,3);

    threepoint(1,2,3) = prob1;
    threepoint(1,3,2) = threepoint(1,2,3);
    threepoint(2,1,3) = threepoint(1,2,3);
    threepoint(2,3,1) = threepoint(1,2,3);
    threepoint(3,1,2) = threepoint(1,2,3);
    threepoint(3,2,1) = threepoint(1,2,3);
    threepoint(1,1,1) = mean(spins(1));
    threepoint(2,2,2) = mean(spins(2));
    threepoint(3,3,3) = mean(spins(3));
    threepoint(1,2,2) = corr(spins(1),spins(2));
    threepoint(2,1,2) = threepoint(1,2,2);
    threepoint(2,2,1) = threepoint(1,2,2);
    threepoint(2,1,1) = threepoint(1,2,2);
    threepoint(1,1,2) = threepoint(1,2,2);
    threepoint(1,2,1) = threepoint(1,2,2);

    threepoint(1,3,3) = corr(spins(1),spins(3));
    threepoint(3,1,3) = threepoint(1,3,3);
    threepoint(3,3,1) = threepoint(1,3,3);
    threepoint(3,1,1) = threepoint(1,3,3);
    threepoint(1,1,3) = threepoint(1,3,3);
    threepoint(1,3,1) = threepoint(1,3,3);

    threepoint(3,2,2) = corr(spins(2),spins(3));
    threepoint(2,3,2) = threepoint(3,2,2);
    threepoint(2,2,3) = threepoint(3,2,2);
    threepoint(2,3,3) = threepoint(3,2,2);
    threepoint(3,3,2) = threepoint(3,2,2);
    threepoint(3,2,3) = threepoint(3,2,2);
end

function [deltaH]  = pointthreeMutual(J, h, mean, corr, numSpins, spins, pair)
    % Calculate the DKL drop for a set of spins connecting to a pair of
    % spins

    deltaH = zeros(numSpins,1);
    
    for i = spins
        j = pair(1);
        k = pair(2);
        if i == j || i == k
            deltaH(i, ind) = 0; 
            continue
        end

        % ----- %
        m = mean([i,j,k]);
        C = corr([i,j,k],[i,j,k]);
        
        step_size = 1/2;
        [h1, J12, J13] = inverse_Ising_GSP_01_helper(m, C, step_size);

        h(i) = h1;
        Jij = J12;
        Jik = J13;


        hjold = h(j);
        hkold = h(k);
        Jjkold = J(j,k);

        hj = h(j) + log(exp(h(i))+1) - log(exp(Jij + h(i)) + 1);
        hk = h(k) + log(exp(h(i))+1) - log(exp(Jik + h(i)) + 1);
        Jjk = J(j,k) - log(exp(h(i)) + 1 ) + log( exp(Jij + h(i))+ 1) + log( exp(Jik + h(i))+ 1) - log( exp(Jij + Jik + h(i))+ 1) ;


        % ----_ %

        deltaH(i) = -log(exp(h(i))+1) + Jij*corr(i,j) + Jik*corr(i,k) + h(i)*mean(i) - mean(i)*log(mean(i))- (1-mean(i))*log(1-mean(i)) + (Jjk - Jjkold)*corr(j,k) + (hj - hjold)*mean(j) + (hk- hkold)*mean(k); 
    end

end

function [h, J,Ent, NewEnt] = find_GSP_update(mean, corr)
    % Updated with new inverse calculator
    % Updated with new change in entropy calculation

    Ent = sum(Entropy(mean)); % Entropy of independent model
    NewEnt = Ent;

    numSpins = length(mean);
    

    h = log(mean./(1-mean));
    J = zeros(numSpins, numSpins);
    binaryJ = zeros(numSpins,numSpins);
    
    MIs_temp = MI2(mean,corr);

    [M,I] = max(MIs_temp, [], 'all', 'linear');
    Ent = Ent - M;
    NewEnt = NewEnt - M;
    [start, ind] = ind2sub(size(corr),I);

    h(ind) = log((mean(ind) - corr(ind, start))./(1 + corr(ind, start) - mean(start) - mean(ind)));
    

    weight = log(corr(ind,start)./(mean(start) - corr(ind,start))) - h(ind);

    h(start) = h(start) + log(exp(h(ind))+ 1) - log(exp(weight + h(ind))+1);


    J(start, ind) = weight;
    J(ind, start) = J(start,ind);
    binaryJ(start, ind) = 1;
    binaryJ(ind,start) = 1;
    

    
    % Initilize list of added pairs
    uniquePairsfull = zeros(2*numSpins - 3,2);

    deltaHs = zeros(numSpins,2*numSpins - 3);

    pair = sort([start, ind]);
    connected_spins = pair;

    spins = 1:numSpins;
    spins(pair) = [];
    [deltaH]  = pointthreeMutual(J, h, mean, corr, numSpins, spins, pair);

    deltaHs(:,1) = deltaH; 
    uniquePairsfull(1,1) = pair(1);
    uniquePairsfull(1,2) = pair(2);
    
    count = 2;
    while count < numSpins

        [M,I] = max(deltaHs, [],"all", "linear");

        dim = size(deltaHs);
        [i,col] = ind2sub(dim,I);


        jk = uniquePairsfull(col,:);
        j = jk(1);
        k = jk(2);

        %-----%

        m = mean([i,j,k]);
        C = corr([i,j,k],[i,j,k]);
        
        step_size = 1;
        [h1, J12, J13] = inverse_Ising_GSP_01_helper(m, C, step_size);

        h(i) = h1;
        Jij = J12;
        Jik = J13;

        %-----%

        J(i, j) = Jij;
        J(j, i) = Jij;
        J(i, k) = Jik;
        J(k, i) = Jik;
        binaryJ(i, j) = 1;
        binaryJ(j, i) = 1;
        binaryJ(i, k) = 1;
        binaryJ(k, i) = 1;


        hjold = h(j);
        hkold = h(k);
        Jjkold = J(j,k);

        hj = h(j) + log(exp(h(i))+1) - log(exp(Jij + h(i)) + 1);
        hk = h(k) + log(exp(h(i))+1) - log(exp(Jik + h(i)) + 1);
        Jjk = J(j,k) - log(exp(h(i)) + 1 ) + log( exp(Jij + h(i))+ 1) + log( exp(Jik + h(i))+ 1) - log( exp(Jij + Jik + h(i))+ 1) ;


        h(j) = hj;
        h(k) = hk;
        J(j,k) = Jjk;
        J(k,j) = J(j,k);


        Ent = Ent - M; % Subtract mutual info from entropy
        deltaH = -log(exp(h(i))+1) + J(i,j)*corr(i,j) + J(i,k)*corr(i,k) + h(i)*mean(i) - mean(i)*log(mean(i))- (1-mean(i))*log(1-mean(i)) + (J(j,k) - Jjkold)*corr(j,k) + (h(j) - hjold)*mean(j) + (h(k)- hkold)*mean(k);

        NewEnt = NewEnt - deltaH;

        
        connected_spins = [connected_spins, i];
    
        spins = 1:numSpins;
        spins(connected_spins) = [];

        % Connect i, j
        pair = sort([i, j]);
        [deltaH]  = pointthreeMutual(J, h, mean, corr, numSpins, spins, pair);

        deltaHs(:,2*(count-1)) = deltaH; 
        uniquePairsfull(2*(count-1),1) = pair(1);
        uniquePairsfull(2*(count-1),2) = pair(2);

        % Connect i, k
        pair = sort([i, k]);
        [deltaH]  = pointthreeMutual(J, h, mean, corr, numSpins, spins, pair);

        deltaHs(:,2*(count-1) + 1) = deltaH; 
        uniquePairsfull(2*(count-1)+1,1) = pair(1);
        uniquePairsfull(2*(count-1)+1,2) = pair(2);

        deltaHs(i, :) = 0;
        count = count + 1
    end

end

function [h, J,Ent, NewEnt] = find_GSP(mean, corr, threepoint)
    % Updated with new inverse calculator

    Ent = sum(Entropy(mean)); % Entropy of independent model
    NewEnt = Ent;

    numSpins = length(mean);


    remaining_spins = 1:numSpins;
    

    h = log(mean./(1-mean));
    J = zeros(numSpins, numSpins);
    binaryJ = zeros(numSpins,numSpins);
    
    MIs_temp = MI2(mean,corr);

    [M,I] = max(MIs_temp, [], 'all', 'linear');
    Ent = Ent - M;
    NewEnt = NewEnt - M;
    [start, ind] = ind2sub(size(corr),I);

    remaining_spins(start) = 0;
    remaining_spins(ind) = 0;

    h(ind) = log((mean(ind) - corr(ind, start))./(1 + corr(ind, start) - mean(start) - mean(ind)));
    

    weight = log(corr(ind,start)./(mean(start) - corr(ind,start))) - h(ind);

    h(start) = h(start) + log(exp(h(ind))+ 1) - log(exp(weight + h(ind))+1);


    J(start, ind) = weight;
    J(ind, start) = J(start,ind);
    binaryJ(start, ind) = 1;
    binaryJ(ind,start) = 1;
    
    
    i = ind;
    j = start;
    if i > j
        j = ind;
        i = start;
    end

    count = 2;
    [MIs, uniquePairs]  = fullthreeMutual(mean, corr, threepoint, numSpins);
    linearind = (i - 1)*numSpins - (i -1)*i/2 + j-i;
    MIs(i, :) = 0;
    MIs(j, :) = 0;

    dim = size(MIs);
    MIstemp = zeros(dim);
    MIstemp(:,linearind) = MIs(:,linearind); % contains all already connected pairs

    while count < numSpins


        %[M,I] = find_max_MIxixjxk_givenMIxixj_xixk(MIs);
        [M,I] = max(MIstemp, [],"all", "linear");

        [i,col] = ind2sub(dim,I);

        MIstemp(i, :) = 0; % remove free node
        MIs(i,:) = 0;

        jk = uniquePairs(col,:);
        j = jk(1);
        k = jk(2);

        %-----%
%         syms u v w
%         Y = vpasolve([mean(i) == (1-mean(j) - mean(k) + corr(j,k))./(1+exp(-u)) + (mean(j) - corr(j,k))./(1+exp(-v-u)) + (mean(k) - corr(j,k))./(1+exp(-w-u)) + (corr(j,k))./(1+exp(-u-v-w)), corr(i,j) == (mean(j) - corr(j,k))./(1+exp(-v-u)) + (corr(j,k))./(1+exp(-u-v-w)), corr(i,k) == (mean(k) - corr(j,k))./(1+exp(-u-w)) + (corr(j,k))./(1+exp(-u-v-w))], [u,v,w]);
% 
%         h(i) = Y.u;
%         Jij = single(Y.v);
%         Jik = single(Y.w);
        m = mean([i,j,k]);
        C = corr([i,j,k],[i,j,k]);
        
        step_size = 1;
        [h1, J12, J13] = inverse_Ising_GSP_01_helper(m, C, step_size);

        h(i) = h1;
        Jij = J12;
        Jik = J13;

        %-----%

        J(i, j) = Jij;
        J(j, i) = Jij;
        J(i, k) = Jik;
        J(k, i) = Jik;
        binaryJ(i, j) = 1;
        binaryJ(j, i) = 1;
        binaryJ(i, k) = 1;
        binaryJ(k, i) = 1;

        remaining_spins(i) = 0;

        hjold = h(j);
        hkold = h(k);
        Jjkold = J(j,k);

        hj = h(j) + log(exp(h(i))+1) - log(exp(Jij + h(i)) + 1);
        hk = h(k) + log(exp(h(i))+1) - log(exp(Jik + h(i)) + 1);
        Jjk = J(j,k) - log(exp(h(i)) + 1 ) + log( exp(Jij + h(i))+ 1) + log( exp(Jik + h(i))+ 1) - log( exp(Jij + Jik + h(i))+ 1) ;


        h(j) = hj;
        h(k) = hk;
        J(j,k) = Jjk;
        J(k,j) = J(j,k);


        xleft = i;
        xright = j;
        if i > j
            xleft = j;
            xright = i;   
        end
        linearind = (xleft - 1)*numSpins - (xleft -1)*xleft/2 + xright - xleft;

        MIstemp(:,linearind) = MIs(:,linearind);

        xleft = i;
        xright = k;
        if i > k
            xleft = k;
            xright = i;   
        end
        linearind = (xleft - 1)*numSpins - (xleft -1)*xleft/2 + xright - xleft;

        MIstemp(:,linearind) = MIs(:,linearind);
        
        Ent = Ent - M; % Subtract mutual info from entropy
        deltaH = -log(exp(h(i))+1) + J(i,j)*corr(i,j) + J(i,k)*corr(i,k) + h(i)*mean(i) - mean(i)*log(mean(i))- (1-mean(i))*log(1-mean(i)) + (J(j,k) - Jjkold)*corr(j,k) + (h(j) - hjold)*mean(j) + (h(k)- hkold)*mean(k);

        NewEnt = NewEnt - deltaH;
        count = count + 1;
    end

end

function [Jguess, hguess, Ent] = FindTree(mean, corr)
    % Find the best tree
    Ent = sum(Entropy(mean)); % Entropy of independent model

    numSpins = length(mean);
    deltaS = MI2(mean, corr);
    deltaS = deltaS - diag(diag(deltaS));
    spins = 1:numSpins;
    Jguess = zeros(numSpins,numSpins);

    hguess = log(mean./(1-mean));
    
    dims = size(corr);

    [M,I] = max(deltaS, [], 'all', 'linear');
        
    [row, col] = ind2sub(dims,I);
    
    %Initialize first connection

    deltaS(row, col) = 0;
    deltaS(col, row) = 0;
    hguess(row) = log( (mean(row) -corr(row,col))/ (1 + corr(row,col) - mean(row) -mean(col)) );
    Jguess(row,col) = log( corr(row,col)/ (mean(col) - corr(row,col))) - hguess(row);
    Jguess(col, row) = Jguess(row, col);

    hguess(col) = hguess(col) + log(exp(hguess(row)) + 1)- log(exp(Jguess(row,col) + hguess(row)  )+ 1);
    numSpins = numSpins -2;
    connectedSpins = [row, col];
    spins = spins(spins~=row );
    spins = spins(spins~=col );
    Ent = Ent - M;
    count = 1
    for i= 1:numSpins
        [M, row, col, trow, tcol] = find_max_in_selected(deltaS, spins, connectedSpins); % Find maximum change in entropy
        Ent = Ent - M;
        connectedSpins = [connectedSpins, row]; %Update list of already connected spins
        spins(trow) = [];
  
        deltaS(row, col) = 0;
        deltaS(col, row) = 0;
        hguess(row) = log( (mean(row) -corr(row,col))/ (1 + corr(row,col) - mean(row) -mean(col)) );
        Jguess(row,col) = log( corr(row,col)/ (mean(col) - corr(row,col))) - hguess(row);
        Jguess(col, row) = Jguess(row, col);

        hguess(col) = hguess(col) + log(exp(hguess(row)) + 1)- log(exp(Jguess(row,col) + hguess(row)  )+ 1);
        count = count + 1
    end
    
end

function [maxValue, rowIndex, colIndex, truerowIndex, truecolIndex] = find_max_in_selected(matrix, selectedRows, selectedCols)
    % Extract the submatrix based on the specified rows and columns
    submatrix = matrix(selectedRows, selectedCols);
    
    % Find the maximum value and its indices in the submatrix
    [maxValue, linearIndex] = max(submatrix(:));
    [rowIndex, colIndex] = ind2sub(size(submatrix), linearIndex);
    
    truerowIndex = rowIndex;
    truecolIndex = colIndex;

    % Map the indices to the original matrix
    rowIndex = selectedRows(rowIndex);
    colIndex = selectedCols(colIndex);
end

function [MIs, uniquePairs]  = fullthreeMutualapprox(mean, corr, numSpins)

    spins = 1:numSpins;

    uniquePairs = zeros(numSpins*(numSpins-1)/2, 2);
    for i = 1:numSpins-1
        for j = i+1:numSpins
            uniquePairs( (i - 1)*numSpins - (i -1)*i/2 + j-i, 1) = i;
            uniquePairs((i - 1)*numSpins - (i -1)*i/2 + j-i, 2) = j;
        end
    end
    TwopointMI = MI2(mean,corr);


    MIs = zeros(numSpins,length(uniquePairs));
    for i = spins
        for ind = 1:length(uniquePairs)
            j = uniquePairs(ind,1);
            k = uniquePairs(ind,2);
            if i == j || i == k
                MIs(i, ind) = 0; 
                continue
            end
            MIs(i, ind) = TwopointMI(i, j) + TwopointMI(i, k);
        end
    end
end

function [MIs, uniquePairs]  = fullthreeMutual(mean, corr, threepoint, numSpins)
    spins = 1:numSpins;

    uniquePairs = zeros(numSpins*(numSpins-1)/2, 2);
    for i = 1:numSpins-1
        for j = i+1:numSpins
            uniquePairs( (i - 1)*numSpins - (i -1)*i/2 + j-i, 1) = i;
            uniquePairs((i - 1)*numSpins - (i -1)*i/2 + j-i, 2) = j;
        end
    end


    MIs = zeros(numSpins,length(uniquePairs));
    for i = spins
        for ind = 1:length(uniquePairs)
            j = uniquePairs(ind,1);
            k = uniquePairs(ind,2);
            if i == j || i == k
                MIs(i, ind) = 0; 
                continue
            end
            MIs(i, ind) = Mutal_info(i, j, k, mean, corr, threepoint);
        end
    end
end

function [H] = Entropy(mean)
    H = -mean.*log(mean)-(1-mean).*log(1-mean);
end

function deltaS = MI2(mean,corr)
    
    a = -mean.*log(mean) - (1-mean).*log(1-mean) - (mean.*log(mean)).' - ((1-mean).*log(1-mean)).';
    b = corr.*log(corr);
    c = (mean - corr).*log(mean - corr) + ((mean - corr).*log(mean - corr)).';
    d = (1 - mean - mean.' +corr).*log((1 - mean - mean.' +corr));
    deltaS = real(a + b + c + d);

    for i = 1:length(mean)
        deltaS(i,i) = 0; %-mean(i)*log(mean(i)) - (1-mean(i))*log(1-mean(i));
    end

    deltaS(isnan(deltaS)) = 0;
end

function [MI] = Mutal_info(xi, xj, xk, mean, corr, threepoint)


    Probs = zeros(2,2,2);

    Probs(2,2,2) = threepoint(xi,xj,xk);
    Probs(1,2,2) = corr(xj,xk) - Probs(2,2,2);
    Probs(2,1,2) = corr(xi,xk) - Probs(2,2,2);
    Probs(1,1,2) = mean(xk) - corr(xi,xk) - corr(xj,xk) + Probs(2,2,2);

    Probs(2,2,1) = corr(xi,xj) - Probs(2,2,2);
    Probs(1,2,1) = mean(xj) - corr(xi,xj) - corr(xj,xk) + Probs(2,2,2);
    Probs(2,1,1) = mean(xi) - corr(xi,xj) - corr(xi,xk) + Probs(2,2,2);
    Probs(1,1,1) = 1 - Probs(2,2,2) - Probs(2,2,1) - Probs(2,1,2) - Probs(1,2,2) - Probs(2,1,1) - Probs(1,2,1) - Probs(1,1,2) ;

    Probsxi = sum(Probs,[2,3]);
    Probsxjxk = sum(Probs,1);

%     sum(-Probsxi.*log(Probsxi), 'all', "omitnan") %Entropy of xi
%     sum(Probsxjxk.*log(Probsxjxk./(sum(Probsxjxk,
%     [1,2]).*sum(Probsxjxk,[1,3]))), 'all', "omitnan") %mi of xj xk,
%     is entropy of xj if xk=xj

    MI = real(sum(Probs.*log(Probs./(Probsxi.*Probsxjxk)), 'all', "omitnan"));

end

function [G, D] = decimate_GSP(G0)
    % Input: nxn unweighted adjacency matrix G0 representing a GSP network.
    %
    % Output: Randomly remove nodes of degree 1 or 2. If we remove a node of
    % degree 2 then we place a new edge between the two neighbors of the
    % removed node. Continue this process until no more nodes of degree 1 or 2
    % exist and return the final adjacency matrix G. Also return the order of
    % decimations in the numDec x 3 matrix D, where D(t,1) is the node removed
    % at step t and D(t,2) and D(t,3) are the two nodes from which the node was
    % removed. We note that D(t,3) = 0 if the removed node had degree 1.
    
    % Size of network:
    N = size(G0,1);
    
    % Initialize things:
    G = double(G0 ~= 0);
    D = zeros(N-1,3);
    
    % Find nodes with degree <= 2:
    degs = sum(G);
    inds = intersect(find(degs <= 2), find(degs >= 1));
    
    % Loop until there are no more feasible nodes:
    counter = 1;
    
    while ~isempty(inds) && sum(sum(G)) ~= 0
        
        Dtemp = zeros(1,3);
        
        % Choose node to remove:
        i = inds(randi(length(inds))); % Remove random node
    %     i = inds(1); % Remove first node in list
    %     i = inds(end); % Remove last node in list
        Dtemp(1) = i;
        
        % Neighbors of i:
        Ni = find(G(i,:));
        Dtemp(2) = Ni(1);
        
        % Remove node from network:
        G(i,Ni) = 0;
        G(Ni,i) = 0;
        
        % If degree of i is 2 then connect two neighbors:
        if length(Ni) == 2
            Dtemp(3) = Ni(2);
            
            G(Ni(1), Ni(2)) = 1;
            G(Ni(2), Ni(1)) = 1;
        end
        
        % Compute new feasible nodes to remove:
        degs = sum(G);
        inds = intersect(find(degs <= 2), find(degs >= 1));
        
        D(counter,:) = Dtemp;
        counter = counter + 1;
        
    end
end

function [Js, connum] =  Jpairs(J, D)

    
    Jcon = J ~= 0;
    Jcon = triu(Jcon);


    Js = zeros(size(Jcon)); % I have no idea why i cant update Jcon
    count = 0;
    for iter = 1:length(D)
        
        i = D(iter,1);
        j = D(iter,2);
        k = D(iter,3);

        count = count + 1;
        Js(i, j) = count;
        if k ~=0 
            count = count + 1;
            Js(i, k) = count; 
        end

    end

    Js = Js + Js.';
    connum = count;

end

function [m, C, X, Z] = correlations_GSP_01(J, h)
    % Updated to be more memory efficient for all derivitives 11/3/23
    % Inputs: nxn matirx J and nx1 col vector h of external fields for a system
    % with 0,1 variables. We require that the interaction matrix J represents a
    % GSP network. We also take as input the order of node decimations (as
    % output by 'decimate_GSP'), where decimations(t,1) is the t^th node to be 
    % decimated and decimations(t,[2,3]) are its two neighbors. We note that
    % decimations(t,3) = 0 if the t^th node only has one neighbor.
    %
    % Output: nx1 vector m of magnetizations, calculated in two steps: (i)
    % decimate nodes of degree 1 or 2 down to a single node, and (ii)
    % recursivley calculate the magnetizations back out in reverse order. We
    % also compute the correlations C(i,j) between all nodes i and j and the
    % partition function Z. See Report 2 for details.
    
    % NOTE: The difference between this function and 'correlations_GSP' is that
    % here consider systems with 0,1 variables.
    
    % Number of nodes:
    n = length(h);

    
    
    %Compute decimation order:
    [Jdec, D] = decimate_GSP(J);
    if sum(sum(Jdec)) ~= 0
        error('Network was not able to be decimated!');
    end
    [Jcon, connum] =  Jpairs(J, D);
    
%     [child, parent1, parent2, h_eff, J_eff] = Decimate(J, h);
% 
% 
%     D = zeros(length(child)-1,3);
% 
%     D(:,1) = child(1:end-1);
%     D(:,2) = parent1(1:end-1);
%     D(:,3) = parent2(1:end-1);

%     c = 0 ;
%     for i = child
%         c = c + log(exp(h_eff(i))+1);
%     end

    % Initialize effective parameters:
    J_eff = J;
    h_eff = h;
    c = 0;
    
    % Loop through node decimations:
    for t = 1:size(D,1)
        
        % Node to decimate and its neighbors:
        i = D(t,1);
        j = D(t,2);
        k = D(t,3);
        
        % Update effective parameters on first neighbor:
        c = c + log(exp(h_eff(i)) + 1);
        h_eff(j) = h_eff(j) - log(exp(h_eff(i)) + 1) + log(exp(J_eff(i,j) + h_eff(i)) + 1);
        
        % If node i has two neighbors:
        if k ~= 0
            
            h_eff(k) = h_eff(k) - log(exp(h_eff(i)) + 1) + log(exp(J_eff(i,k) + h_eff(i)) + 1);
            J_eff(j,k) = J_eff(j,k) + log(exp(h_eff(i)) + 1) - log(exp(J_eff(i,j) + h_eff(i)) + 1)...
                - log(exp(J_eff(i,k) + h_eff(i)) + 1) + log(exp(J_eff(i,j) + J_eff(i,k) + h_eff(i)) + 1);
            J_eff(k,j) = J_eff(j,k);
            
        end
    end
    
    % Things we are going to need to compute:
    m = zeros(n,1); % Magnetizations
    C = zeros(n); % Correlations between nodes that interact
    dm_dh = zeros(n); % Susceptibilities
    dC_dh = zeros(connum,n); % Derivatives of correlations with respect to external fields
    dh_dh = eye(n); % Derivatives of external fields with respect to external fields
    dJ_dh = zeros(connum,n); % Derivatives of interactions with respect to external fields
    dh_dJ = zeros(n,connum); % Derivatives of external fields with respect to interactions
    Jpair = zeros(connum);
    %dJ_dJ = zeros(n,n,n,n); % Derivatives of interactions with respect to interactions
    
    % Make self-derivatives unity:
    for ind1 = 1:n
        for ind2 = 1:n
            if Jcon(ind1,ind2) ~= 0
                Jpair(Jcon(ind1,ind2),Jcon(ind1,ind2)) = 1;
                Jpair(Jcon(ind2,ind1),Jcon(ind1,ind2)) = 1;
                Jpair(Jcon(ind1,ind2),Jcon(ind2,ind1)) = 1;
                Jpair(Jcon(ind2,ind1),Jcon(ind2,ind1)) = 1;
            end
%             dJ_dJ(ind1,ind2,ind1,ind2) = 1;
%             dJ_dJ(ind1,ind2,ind2,ind1) = 1;
%             dJ_dJ(ind2,ind1,ind1,ind2) = 1;
%             dJ_dJ(ind2,ind1,ind2,ind1) = 1;
        end
    end
    
    % Compute things for final node:
    i0 = j ;%child(end);
    m(i0) = 1/(1 + exp(-h_eff(i0)));
    C(i0,i0) = m(i0);
    dm_dh(i0,i0) = exp(-h_eff(i0))/(1 + exp(-h_eff(i0)))^2;
    
    % Compute partition function:
    Z = exp(c)*(exp(h_eff(i0)) + 1);
    
    % Loop over nodes in reverse order from which they were decimated:
    for t_j = size(D,1):-1:1
        
        % Node to take derivative with respect to and its decimation neighbors:
        j1 = D(t_j,1);
        j2 = D(t_j,2);
        j3 = D(t_j,3);
        
        % Compute magnetizations and derivatives:
        
        % If j1 only has one neighbor:
        if j3 == 0
            
            % Compute magnetizations and correlations:
            m(j1) = (1 - m(j2))/(1 + exp(-h_eff(j1))) + m(j2)/(1 + exp(-J_eff(j1,j2) - h_eff(j1)));
            C(j1,j1) = m(j1);
            
            C(j1,j2) = m(j2)/(1 + exp(-J_eff(j1,j2) - h_eff(j1)));
            C(j2,j1) = C(j1,j2);
            
            % Compute derivatives of h(j2) with respect to h(j1) and J(j1,j2):
            dh_dh(j2,j1) = -1/(1 + exp(-h_eff(j1))) + 1/(1 + exp(-J_eff(j1,j2) - h_eff(j1)));
            dh_dJ(j2,Jcon(j1,j2)) = 1/(1 + exp(-J_eff(j1,j2) - h_eff(j1)));
            dh_dJ(j2,Jcon(j2,j1)) = dh_dJ(j2,Jcon(j1,j2));
            
            % Compute derivative of external field of final node with respect
            % to h(j1) and J(j1,j2). Dependence goes through h(j2):
            dh_dh(i0,j1) = dh_dh(i0,j2)*dh_dh(j2,j1);
            dh_dJ(i0,Jcon(j1,j2)) = dh_dh(i0,j2)*dh_dJ(j2,Jcon(j1,j2));
            dh_dJ(i0,Jcon(j2,j1)) = dh_dJ(i0,Jcon(j1,j2));
            
            % Loop over nodes between final node and j1:
            for t_i = size(D,1):-1:(t_j + 1)
                
                % Node to take derivative of and its decimation neighbors:
                i1 = D(t_i,1);
                i2 = D(t_i,2);
                i3 = D(t_i,3);
                
                % Compute derivatives of h(i1), J(i1,i2), and J(i1,i3) with
                % respect to h(j1) and J(j1,j2). All dependencies go through h(j2):
                dh_dh(i1,j1) = dh_dh(i1,j2)*dh_dh(j2,j1);
                dh_dJ(i1,Jcon(j1,j2)) = dh_dh(i1,j2)*dh_dJ(j2,Jcon(j1,j2));
                dh_dJ(i1,Jcon(j2,j1)) = dh_dJ(i1,Jcon(j1,j2));
                
                dJ_dh(Jcon(i1,i2),j1) = dJ_dh(Jcon(i1,i2),j2)*dh_dh(j2,j1);
                dJ_dh(Jcon(i2,i1),j1) = dJ_dh(Jcon(i1,i2),j1);


                Jpair(Jcon(i1,i2),Jcon(j1,j2)) = dJ_dh(Jcon(i1,i2),j2)*dh_dJ(j2,Jcon(j1,j2));
                Jpair(Jcon(i2,i1),Jcon(j1,j2)) = Jpair(Jcon(i1,i2),Jcon(j1,j2));
                Jpair(Jcon(i1,i2),Jcon(j2,j1)) = Jpair(Jcon(i1,i2),Jcon(j1,j2));
                Jpair(Jcon(i2,i1),Jcon(j2,j1)) = Jpair(Jcon(i1,i2),Jcon(j1,j2));

%                 dJ_dJ(i1,i2,j1,j2) = dJ_dh(i1,i2,j2)*dh_dJ(j2,j1,j2);
%                 dJ_dJ(i2,i1,j1,j2) = dJ_dJ(i1,i2,j1,j2);
%                 dJ_dJ(i1,i2,j2,j1) = dJ_dJ(i1,i2,j1,j2);
%                 dJ_dJ(i2,i1,j2,j1) = dJ_dJ(i1,i2,j1,j2);
%                 
                % If i1 has two neighbors:
                if i3 ~= 0
                    
                    dJ_dh(Jcon(i1,i3),j1) = dJ_dh(Jcon(i1,i3),j2)*dh_dh(j2,j1);
                    dJ_dh(Jcon(i3,i1),j1) = dJ_dh(Jcon(i1,i3),j1);

                    Jpair(Jcon(i1,i3),Jcon(j1,j2)) = dJ_dh(Jcon(i1,i3),j2)*dh_dJ(j2,Jcon(j1,j2));
                    Jpair(Jcon(i3,i1),Jcon(j1,j2)) = Jpair(Jcon(i1,i3),Jcon(j1,j2));
                    Jpair(Jcon(i1,i3),Jcon(j2,j1)) = Jpair(Jcon(i1,i3),Jcon(j1,j2));
                    Jpair(Jcon(i3,i1),Jcon(j2,j1)) = Jpair(Jcon(i1,i3),Jcon(j1,j2));

%                     dJ_dJ(i1,i3,j1,j2) = dJ_dh(i1,i3,j2)*dh_dJ(j2,j1,j2);
%                     dJ_dJ(i3,i1,j1,j2) = dJ_dJ(i1,i3,j1,j2);
%                     dJ_dJ(i1,i3,j2,j1) = dJ_dJ(i1,i3,j1,j2);
%                     dJ_dJ(i3,i1,j2,j1) = dJ_dJ(i1,i3,j1,j2);
                    
                end
      
            end
            
        % If j1 has two neighbors:
        else
            
            % Compute magnetizations and correlations:
            m(j1) = (1 - m(j2) - m(j3) + C(j2,j3))/(1 + exp(-h_eff(j1)))...
                + (m(j2) - C(j2,j3))/(1 + exp(-J_eff(j1,j2) - h_eff(j1)))...
                + (m(j3) - C(j2,j3))/(1 + exp(-J_eff(j1,j3) - h_eff(j1)))...
                + C(j2,j3)/(1 + exp(-J_eff(j1,j2) - J_eff(j1,j3) - h_eff(j1)));
            C(j1,j1) = m(j1);
            
            C(j1,j2) = (m(j2) - C(j2,j3))/(1 + exp(-J_eff(j1,j2) - h_eff(j1)))...
                + C(j2,j3)/(1 + exp(-J_eff(j1,j2) - J_eff(j1,j3) - h_eff(j1)));
            C(j2,j1) = C(j1,j2);
            
            C(j1,j3) = (m(j3) - C(j2,j3))/(1 + exp(-J_eff(j1,j3) - h_eff(j1)))...
                + C(j2,j3)/(1 + exp(-J_eff(j1,j2) - J_eff(j1,j3) - h_eff(j1)));
            C(j3,j1) = C(j1,j3);
            
            % Compute derivatives of h(j2), h(j3), and J(j2,j3) with respect to
            % h(j1), J(j1,j2), and J(j1,j3):
            dh_dh(j2,j1) = -1/(1 + exp(-h_eff(j1))) + 1/(1 + exp(-J_eff(j1,j2) - h_eff(j1)));
            dh_dJ(j2,Jcon(j1,j2)) = 1/(1 + exp(-J_eff(j1,j2) - h_eff(j1)));
            dh_dJ(j2,Jcon(j2,j1)) = dh_dJ(j2,Jcon(j1,j2));
            
            dh_dh(j3,j1) = -1/(1 + exp(-h_eff(j1))) + 1/(1 + exp(-J_eff(j1,j3) - h_eff(j1)));
            dh_dJ(j3,Jcon(j1,j3)) = 1/(1 + exp(-J_eff(j1,j3) - h_eff(j1)));
            dh_dJ(j3,Jcon(j3,j1)) = dh_dJ(j3,Jcon(j1,j3));
            
            dJ_dh(Jcon(j2,j3),j1) = 1/(1 + exp(-h_eff(j1))) - 1/(1 + exp(-J_eff(j1,j2) - h_eff(j1)))...
                - 1/(1 + exp(-J_eff(j1,j3) - h_eff(j1))) + 1/(1 + exp(-J_eff(j1,j2) - J_eff(j1,j3) - h_eff(j1)));
            dJ_dh(Jcon(j3,j2),j1) = dJ_dh(Jcon(j2,j3),j1);


            Jpair(Jcon(j2,j3),Jcon(j1,j2)) = -1/(1 + exp(-J_eff(j1,j2) - h_eff(j1)))...
                + 1/(1 + exp(-J_eff(j1,j2) - J_eff(j1,j3) - h_eff(j1)));
            Jpair(Jcon(j2,j3),Jcon(j2,j1)) = Jpair(Jcon(j2,j3),Jcon(j1,j2));
            Jpair(Jcon(j3,j2),Jcon(j1,j2)) = Jpair(Jcon(j2,j3),Jcon(j1,j2));
            Jpair(Jcon(j3,j2),Jcon(j2,j1)) = Jpair(Jcon(j2,j3),Jcon(j1,j2));


%             dJ_dJ(j2,j3,j1,j2) = -1/(1 + exp(-J_eff(j1,j2) - h_eff(j1)))...
%                 + 1/(1 + exp(-J_eff(j1,j2) - J_eff(j1,j3) - h_eff(j1)));
%             dJ_dJ(j2,j3,j2,j1) = dJ_dJ(j2,j3,j1,j2);
%             dJ_dJ(j3,j2,j1,j2) = dJ_dJ(j2,j3,j1,j2);
%             dJ_dJ(j3,j2,j2,j1) = dJ_dJ(j2,j3,j1,j2);


            Jpair(Jcon(j2,j3),Jcon(j1,j3)) = -1/(1 + exp(-J_eff(j1,j3) - h_eff(j1)))...
                + 1/(1 + exp(-J_eff(j1,j2) - J_eff(j1,j3) - h_eff(j1)));
            Jpair(Jcon(j2,j3),Jcon(j3,j1)) = Jpair(Jcon(j2,j3),Jcon(j1,j3));
            Jpair(Jcon(j3,j2),Jcon(j1,j3)) = Jpair(Jcon(j2,j3),Jcon(j1,j3));
            Jpair(Jcon(j3,j2),Jcon(j3,j1)) = Jpair(Jcon(j2,j3),Jcon(j1,j3));

%             dJ_dJ(j2,j3,j1,j3) = -1/(1 + exp(-J_eff(j1,j3) - h_eff(j1)))...
%                 + 1/(1 + exp(-J_eff(j1,j2) - J_eff(j1,j3) - h_eff(j1)));
%             dJ_dJ(j2,j3,j3,j1) = dJ_dJ(j2,j3,j1,j3);
%             dJ_dJ(j3,j2,j1,j3) = dJ_dJ(j2,j3,j1,j3);
%             dJ_dJ(j3,j2,j3,j1) = dJ_dJ(j2,j3,j1,j3);
            
            % Compute derivative of external field of final node with respect
            % to h(j1), J(j1,j2), and J(j1,j3). Dependencies go through h(j2),
            % h(j3), and J(j2,j3):
            dh_dh(i0,j1) = dh_dh(i0,j2)*dh_dh(j2,j1) + dh_dh(i0,j3)*dh_dh(j3,j1)...
                + dh_dJ(i0,Jcon(j2,j3))*dJ_dh(Jcon(j2,j3),j1);
            dh_dJ(i0,Jcon(j1,j2)) = dh_dh(i0,j2)*dh_dJ(j2,Jcon(j1,j2)) + dh_dJ(i0,Jcon(j2,j3))*Jpair(Jcon(j2,j3),Jcon(j1,j2));  %dJ_dJ(j2,j3,j1,j2);
            dh_dJ(i0,Jcon(j2,j1)) = dh_dJ(i0,Jcon(j1,j2));
            dh_dJ(i0,Jcon(j1,j3)) = dh_dh(i0,j3)*dh_dJ(j3,Jcon(j1,j3)) + dh_dJ(i0,Jcon(j2,j3))*Jpair(Jcon(j2,j3),Jcon(j1,j3));  %dJ_dJ(j2,j3,j1,j3);
            dh_dJ(i0,Jcon(j3,j1)) = dh_dJ(i0,Jcon(j1,j3));
            
            % Loop over nodes between final node and j1:
            for t_i = size(D,1):-1:(t_j + 1)
                
                % Node to take derivative of and its decimation neighbors:
                i1 = D(t_i,1);
                i2 = D(t_i,2);
                i3 = D(t_i,3);
                
                % Compute derivatives of h(i1), J(i1,i2), and J(i1,i3) with
                % respect to h(j1), J(j1,j2) and J(j1,j3). All dependencies go
                % through h(j2), h(j3), and J(j2,j3):
                dh_dh(i1,j1) = dh_dh(i1,j2)*dh_dh(j2,j1) + dh_dh(i1,j3)*dh_dh(j3,j1)...
                    + dh_dJ(i1,Jcon(j2,j3))*dJ_dh(Jcon(j2,j3),j1);
                dh_dJ(i1,Jcon(j1,j2)) = dh_dh(i1,j2)*dh_dJ(j2,Jcon(j1,j2)) + dh_dJ(i1,Jcon(j2,j3))*Jpair(Jcon(j2,j3),Jcon(j1,j2)); %dJ_dJ(j2,j3,j1,j2);
                dh_dJ(i1,Jcon(j2,j1)) = dh_dJ(i1,Jcon(j1,j2));
                dh_dJ(i1,Jcon(j1,j3)) = dh_dh(i1,j3)*dh_dJ(j3,Jcon(j1,j3)) + dh_dJ(i1,Jcon(j2,j3))*Jpair(Jcon(j2,j3),Jcon(j1,j3)); %dJ_dJ(j2,j3,j1,j3);
                dh_dJ(i1,Jcon(j3,j1)) = dh_dJ(i1,Jcon(j1,j3));
                
                dJ_dh(Jcon(i1,i2),j1) = dJ_dh(Jcon(i1,i2),j2)*dh_dh(j2,j1) + dJ_dh(Jcon(i1,i2),j3)*dh_dh(j3,j1)...
                    + dJ_dh(Jcon(j2,j3),j1)*Jpair(Jcon(i1,i2),Jcon(j2,j3)); %dJ_dJ(i1,i2,j2,j3)*
                dJ_dh(Jcon(i2,i1),j1) = dJ_dh(Jcon(i1,i2),j1);

                Jpair(Jcon(i1,i2),Jcon(j1,j2)) = dJ_dh(Jcon(i1,i2),j2)*dh_dJ(j2,Jcon(j1,j2)) + Jpair(Jcon(i1,i2),Jcon(j2,j3))*Jpair(Jcon(j2,j3),Jcon(j1,j2));
                Jpair(Jcon(i1,i2),Jcon(j2,j1)) = Jpair(Jcon(i1,i2),Jcon(j1,j2));
                Jpair(Jcon(i2,i1),Jcon(j1,j2)) = Jpair(Jcon(i1,i2),Jcon(j1,j2));
                Jpair(Jcon(i2,i1),Jcon(j2,j1)) = Jpair(Jcon(i1,i2),Jcon(j1,j2));

%                 dJ_dJ(i1,i2,j1,j2) = dJ_dh(i1,i2,j2)*dh_dJ(j2,j1,j2) + dJ_dJ(i1,i2,j2,j3)*dJ_dJ(j2,j3,j1,j2);
%                 dJ_dJ(i1,i2,j2,j1) = dJ_dJ(i1,i2,j1,j2);
%                 dJ_dJ(i2,i1,j1,j2) = dJ_dJ(i1,i2,j1,j2);
%                 dJ_dJ(i2,i1,j2,j1) = dJ_dJ(i1,i2,j1,j2);

                Jpair(Jcon(i1,i2),Jcon(j1,j3)) = dJ_dh(Jcon(i1,i2),j3)*dh_dJ(j3,Jcon(j1,j3)) + Jpair(Jcon(i1,i2),Jcon(j2,j3))*Jpair(Jcon(j2,j3),Jcon(j1,j3));
                Jpair(Jcon(i1,i2),Jcon(j3,j1)) = Jpair(Jcon(i1,i2),Jcon(j1,j3));
                Jpair(Jcon(i2,i1),Jcon(j1,j3)) = Jpair(Jcon(i1,i2),Jcon(j1,j3));
                Jpair(Jcon(i2,i1),Jcon(j3,j1)) = Jpair(Jcon(i1,i2),Jcon(j1,j3));

%                 dJ_dJ(i1,i2,j1,j3) = dJ_dh(i1,i2,j3)*dh_dJ(j3,j1,j3) + dJ_dJ(i1,i2,j2,j3)*dJ_dJ(j2,j3,j1,j3);
%                 dJ_dJ(i1,i2,j3,j1) = dJ_dJ(i1,i2,j1,j3);
%                 dJ_dJ(i2,i1,j1,j3) = dJ_dJ(i1,i2,j1,j3);
%                 dJ_dJ(i2,i1,j3,j1) = dJ_dJ(i1,i2,j1,j3);
                
                % If i1 has two neighbors:
                if i3 ~= 0
                    
                    dJ_dh(Jcon(i1,i3),j1) = dJ_dh(Jcon(i1,i3),j2)*dh_dh(j2,j1) + dJ_dh(Jcon(i1,i3),j3)*dh_dh(j3,j1)...
                        + dJ_dh(Jcon(j2,j3),j1)*Jpair(Jcon(i1,i3),Jcon(j2,j3)); %dJ_dJ(i1,i3,j2,j3)*
                    dJ_dh(Jcon(i3,i1),j1) = dJ_dh(Jcon(i1,i3),j1);

                    Jpair(Jcon(i1,i3),Jcon(j1,j2)) = dJ_dh(Jcon(i1,i3),j2)*dh_dJ(j2,Jcon(j1,j2)) + Jpair(Jcon(i1,i3),Jcon(j2,j3))*Jpair(Jcon(j2,j3),Jcon(j1,j2));
                    Jpair(Jcon(i1,i3),Jcon(j2,j1)) = Jpair(Jcon(i1,i3),Jcon(j1,j2));
                    Jpair(Jcon(i3,i1),Jcon(j1,j2)) = Jpair(Jcon(i1,i3),Jcon(j1,j2));
                    Jpair(Jcon(i3,i1),Jcon(j2,j1)) = Jpair(Jcon(i1,i3),Jcon(j1,j2));

%                     dJ_dJ(i1,i3,j1,j2) = dJ_dh(i1,i3,j2)*dh_dJ(j2,j1,j2) + dJ_dJ(i1,i3,j2,j3)*dJ_dJ(j2,j3,j1,j2);
%                     dJ_dJ(i1,i3,j2,j1) = dJ_dJ(i1,i3,j1,j2);
%                     dJ_dJ(i3,i1,j1,j2) = dJ_dJ(i1,i3,j1,j2);
%                     dJ_dJ(i3,i1,j2,j1) = dJ_dJ(i1,i3,j1,j2);


                    Jpair(Jcon(i1,i3),Jcon(j1,j3)) = dJ_dh(Jcon(i1,i3),j3)*dh_dJ(j3,Jcon(j1,j3)) + Jpair(Jcon(i1,i3),Jcon(j2,j3))*Jpair(Jcon(j2,j3),Jcon(j1,j3));
                    Jpair(Jcon(i1,i3),Jcon(j3,j1)) = Jpair(Jcon(i1,i3),Jcon(j1,j3));
                    Jpair(Jcon(i3,i1),Jcon(j1,j3)) = Jpair(Jcon(i1,i3),Jcon(j1,j3));
                    Jpair(Jcon(i3,i1),Jcon(j3,j1)) = Jpair(Jcon(i1,i3),Jcon(j1,j3));
% 
%                     dJ_dJ(i1,i3,j1,j3) = dJ_dh(i1,i3,j3)*dh_dJ(j3,j1,j3) + dJ_dJ(i1,i3,j2,j3)*dJ_dJ(j2,j3,j1,j3);
%                     dJ_dJ(i1,i3,j3,j1) = dJ_dJ(i1,i3,j1,j3);
%                     dJ_dJ(i3,i1,j1,j3) = dJ_dJ(i1,i3,j1,j3);
%                     dJ_dJ(i3,i1,j3,j1) = dJ_dJ(i1,i3,j1,j3);
                
                end
                
            end
            
        end
        
        % Compute susceptibilities:
        
        % Compute susceptibility with final node (derivative of m(i0) with
        % respect to h(j1)). Dependence goes through h(i0):
        dm_dh(i0,j1) = dm_dh(i0,i0)*dh_dh(i0,j1);
        dm_dh(j1,i0) = dm_dh(i0,j1);
        
        % Loop over nodes between final node and j1:
        for t_i = size(D,1):-1:(t_j + 1)
            
            % Node to take derivative of and its decimation neighbors:
            i1 = D(t_i,1);
            i2 = D(t_i,2);
            i3 = D(t_i,3);
                
            % If i1 only has one neighbor:
            if i3 == 0
                
                % Useful quantities:
                l_h = 1/(1 + exp(-h_eff(i1)));
                l_hJ = 1/(1 + exp(-J_eff(i1,i2) - h_eff(i1)));
                dl_h = exp(-h_eff(i1))/(1 + exp(-h_eff(i1)))^2;
                dl_hJ = exp(-J_eff(i1,i2) - h_eff(i1))/(1 + exp(-J_eff(i1,i2) - h_eff(i1)))^2;
                
                % Compute susceptibility (derivative of m(i1) with respect to h(j1)):
                dm_dh(i1,j1) = dl_h*(1 - m(i2))*dh_dh(i1,j1) + dl_hJ*m(i2)*(dh_dh(i1,j1) + dJ_dh(Jcon(i1,i2),j1))...
                    + (-l_h + l_hJ)*dm_dh(i2,j1);
                dm_dh(j1,i1) = dm_dh(i1,j1);
                
                % Compute derivative of C(i1,i2) with respect to h(j1):
                dC_dh(Jcon(i1,i2),j1) = dl_hJ*m(i2)*(dh_dh(i1,j1) + dJ_dh(Jcon(i1,i2),j1)) + l_hJ*dm_dh(i2,j1);
                dC_dh(Jcon(i2,i1),j1) = dC_dh(Jcon(i1,i2),j1);
                
            % If i1 has two neighbors:
            else
                
                % Useful quantities:
                l_h = 1/(1 + exp(-h_eff(i1)));
                l_hJ2 = 1/(1 + exp(-J_eff(i1,i2) - h_eff(i1)));
                l_hJ3 = 1/(1 + exp(-J_eff(i1,i3) - h_eff(i1)));
                l_hJJ = 1/(1 + exp(-J_eff(i1,i2) - J_eff(i1,i3) - h_eff(i1)));
                dl_h = exp(-h_eff(i1))/(1 + exp(-h_eff(i1)))^2;
                dl_hJ2 = exp(-J_eff(i1,i2) - h_eff(i1))/(1 + exp(-J_eff(i1,i2) - h_eff(i1)))^2;
                dl_hJ3 = exp(-J_eff(i1,i3) - h_eff(i1))/(1 + exp(-J_eff(i1,i3) - h_eff(i1)))^2;
                dl_hJJ = exp(-J_eff(i1,i2) - J_eff(i1,i3) - h_eff(i1))/(1 + exp(-J_eff(i1,i2) - J_eff(i1,i3) - h_eff(i1)))^2;
                
                % Compute susceptibility (derivative of m(i1) with respect to h(j1)):
                dm_dh(i1,j1) = dl_h*(1 - m(i2) - m(i3) + C(i2,i3))*dh_dh(i1,j1)...
                    + dl_hJ2*(m(i2) - C(i2,i3))*(dh_dh(i1,j1) + dJ_dh(Jcon(i1,i2),j1))...
                    + dl_hJ3*(m(i3) - C(i2,i3))*(dh_dh(i1,j1) + dJ_dh(Jcon(i1,i3),j1))...
                    + dl_hJJ*C(i2,i3)*(dh_dh(i1,j1) + dJ_dh(Jcon(i1,i2),j1) + dJ_dh(Jcon(i1,i3),j1))...
                    + (-l_h + l_hJ2)*dm_dh(i2,j1) + (-l_h + l_hJ3)*dm_dh(i3,j1)...
                    + (l_h - l_hJ2 - l_hJ3 + l_hJJ)*dC_dh(Jcon(i2,i3),j1);
                dm_dh(j1,i1) = dm_dh(i1,j1);
                
                % Compute derivative of C(i1,i2) with respect to h(j1):
                dC_dh(Jcon(i1,i2),j1) = dl_hJ2*(m(i2) - C(i2,i3))*(dh_dh(i1,j1) + dJ_dh(Jcon(i1,i2),j1))...
                    + dl_hJJ*C(i2,i3)*(dh_dh(i1,j1) + dJ_dh(Jcon(i1,i2),j1) + dJ_dh(Jcon(i1,i3),j1))...
                    + l_hJ2*dm_dh(i2,j1) + (-l_hJ2 + l_hJJ)*dC_dh(Jcon(i2,i3),j1);
                dC_dh(Jcon(i2,i1),j1) = dC_dh(Jcon(i1,i2),j1);
                
                % Compute derivative of C(i1,i3) with respect to h(j1):
                dC_dh(Jcon(i1,i3),j1) = dl_hJ3*(m(i3) - C(i2,i3))*(dh_dh(i1,j1) + dJ_dh(Jcon(i1,i3),j1))...
                    + dl_hJJ*C(i2,i3)*(dh_dh(i1,j1) + dJ_dh(Jcon(i1,i2),j1) + dJ_dh(Jcon(i1,i3),j1))...
                    + l_hJ3*dm_dh(i3,j1) + (-l_hJ3 + l_hJJ)*dC_dh(Jcon(i2,i3),j1);
                dC_dh(Jcon(i3,i1),j1) = dC_dh(Jcon(i1,i3),j1);
                
            end
            
        end
        
    end
    
    % Compute correlations between all nodes:
    C = dm_dh + m*m';
    C(logical(eye(n))) = m;
    
    % Susceptibility:
    X = C - m*m';
end

function [mean,corr] = findStats(J,h)
    % Find the statistics of the GSP network defined in the Ising
    % interaction matrix J and external field vector h. Only finds
    % some of the correlation between nodes. Use correlations_GSP_01
    % instead for full correlation matrix. 


    [child, parent1, parent2, heff, Jeff] = Decimate(J, h);

    child = flip(child,2);
    parent1 = flip(parent1,2);
    parent2 = flip(parent2,2); 

    c = 0 ;
    for i = child
        c = c + log(exp(heff(i))+1);
    end

    mean = 1./(1+exp(-heff));

    corr = zeros(size(J));

    dim = size(child, 2);

    xi = child(2);
    xj = parent1(2);
    mean(xi) = (1-mean(xj))./(1+exp(-heff(xi))) + mean(xj)./(1+exp(-Jeff(xi, xj)-heff(xi)));
    corr(xi,xj) = mean(xj)./(1+exp(-Jeff(xi, xj)-heff(xi)));
    corr(xj,xi) = corr(xi,xj);

    for i = 3:dim
    
        xi = child(i);
        xj = parent1(i);
        xk = parent2(i);
        if xk == 0
            mean(xi) = (1-mean(xj))./(1+exp(-heff(xi))) + mean(xj)./(1+exp(-Jeff(xi, xj)-heff(xi)));
            corr(xi,xj) = mean(xj)./(1+exp(-Jeff(xi, xj)-heff(xi)));
            corr(xj,xi) = corr(xi,xj);
        else
            mean(xi) = (1-mean(xj) - mean(xk) + corr(xj,xk))./(1+exp(-heff(xi))) + (mean(xj) - corr(xj,xk))./(1+exp(-Jeff(xi, xj)-heff(xi))) + (mean(xk) - corr(xj,xk))./(1+exp(-Jeff(xi, xk)-heff(xi))) + (corr(xj,xk))./(1+exp(-Jeff(xi, xj)-Jeff(xi, xk)-heff(xi)));
            corr(xi,xj) = (mean(xj) - corr(xj,xk))./(1+exp(-Jeff(xi, xj)-heff(xi))) + (corr(xj,xk))./(1+exp(-Jeff(xi, xj)-Jeff(xi, xk)-heff(xi)));
            corr(xj,xi) = corr(xi,xj);
            corr(xi,xk) = (mean(xk) - corr(xj,xk))./(1+exp(-Jeff(xi, xk)-heff(xi))) + (corr(xj,xk))./(1+exp(-Jeff(xi, xj)-Jeff(xi, xk)-heff(xi)));
            corr(xk,xi) = corr(xi,xk);
        end

    end


end

function [child, parent1, parent2, heff, Jeff] = Decimate(J, h)
    % Decimate the GSP network defined in the connection matrix J
    % Not used in correlations_GSP_01.
    numSpins = length(h);
    G = graph(J, 'upper');
    nodeNames = cellstr(num2str((1:numSpins).')); % Convert indices to cell array of strings
    G.Nodes.Name = nodeNames; % Assign node names
    Jeff = J;

    parent1 = [];
    parent2 = [];
    child = [];
    nodenum = height(G.Nodes);

    heff = h;
    while nodenum > 1
        [m, Gind] = min(G.degree);
        Gneigh = neighbors(G,Gind);
        parents = outedges(G,Gind);
        temp = G.Nodes.Name(Gind);
        ind = str2num(temp{1}); % set as labled node
        child = [child, ind];

        temp = G.Nodes.Name(Gneigh(1));
        neigh(1) = str2num(temp{1}); % set as labled neighbor
        heff(neigh(1)) = heff(neigh(1)) - log(exp(heff(ind)) + 1) + log(exp(heff(ind) + G.Edges.Weight(parents(1))) + 1);
        
        parent1 = [parent1, neigh(1)];
        if m ==2
            temp = G.Nodes.Name(Gneigh(2));
            neigh(2) = str2num(temp{1}); % set as labled neighbor
            heff(neigh(2)) = heff(neigh(2)) - log(exp(heff(ind)) + 1) + log(exp(heff(ind) + G.Edges.Weight(parents(2))) + 1);
            parent2 = [parent2, neigh(2)];
            edgeind = findedge(G,Gneigh(1),Gneigh(2));
            if edgeind == 0
                G = addedge(G,Gneigh(1),Gneigh(2), log(exp(heff(ind)) + 1 ) - log( exp(G.Edges.Weight(parents(1)) + heff(ind))+ 1)- log( exp(G.Edges.Weight(parents(2)) + heff(ind))+ 1) + log( exp(G.Edges.Weight(parents(1))+G.Edges.Weight(parents(2)) + heff(ind))+ 1) ); 
                Jeff(neigh(1), neigh(2)) = log(exp(heff(ind)) + 1 ) - log( exp(G.Edges.Weight(parents(1)) + heff(ind))+ 1)- log( exp(G.Edges.Weight(parents(2)) + heff(ind))+ 1) + log( exp(G.Edges.Weight(parents(1))+G.Edges.Weight(parents(2)) + heff(ind))+ 1);
                Jeff(neigh(2), neigh(1)) = Jeff(neigh(1), neigh(2));
            else
                weight = G.Edges.Weight(edgeind);
                Jeff(neigh(1), neigh(2)) = weight + log(exp(heff(ind)) + 1 ) - log( exp(G.Edges.Weight(parents(1)) + heff(ind))+ 1)- log( exp(G.Edges.Weight(parents(2)) + heff(ind))+ 1) + log( exp(G.Edges.Weight(parents(1))+G.Edges.Weight(parents(2)) + heff(ind))+ 1);
                Jeff(neigh(2), neigh(1)) = Jeff(neigh(1), neigh(2));
                G.Edges.Weight(edgeind) = weight + log(exp(heff(ind)) + 1 ) - log( exp(G.Edges.Weight(parents(1)) + heff(ind))+ 1)- log( exp(G.Edges.Weight(parents(2)) + heff(ind))+ 1) + log( exp(G.Edges.Weight(parents(1))+G.Edges.Weight(parents(2)) + heff(ind))+ 1) ;     
            end
        else
            parent2 = [parent2, 0];
        end

        %remove node
        G = rmnode(G, Gind);
        nodenum = nodenum - 1;

    end

    [m, Gind] = min(G.degree);
    temp = G.Nodes.Name(Gind);
    ind = str2num(temp{1}); % set as labled node
    child = [child, ind];
    parent1 = [parent1,0];
    parent2 = [parent2,0];

end

function [treeleaves, treeconnection, J] = PruneGSP(J)
    % Remove Leaves from GSP network defined in Ising interaction matrix J.
    % Returns the matrix without leaves.
    treeleaves = [];
    treeconnection = [];
    dim = size(J, 1);
    [leaves, nodes] = Findleaves(J);
    while isempty(leaves) == 0
        
        treeleaves = [treeleaves, leaves];
        treeconnection = [treeconnection, nodes];
        J(leaves, : ) = zeros(size(leaves,2), dim);
        J(:, leaves) = zeros(dim, size(leaves,2));
        [leaves, nodes] = Findleaves(J);
    end

    if length(treeleaves) > 1
        if treeconnection(end-1:end) == flip(treeleaves(end-1:end))
    
            treeleaves(end) = [];
            treeconnection(end) = [];
        end
    end
end

function [leaves, nodes] = Findleaves(J)
    % Function takes in a Ising interaction matrix J and determines the
    % nodes that have one parent. Retuns the child nodes in leaves and the
    % parent nodes as nodes. 

    dim = size(J,1); 
    leaves = [];
    nodes = [];
    for i = 1:dim
        if nnz(J(:,i)) == 1
            leaves = [leaves, i];
            nodes = [nodes, find(J(:,i))];
        end
    end

end 

function [J, h] = RandomGSPIsing(numSpins)
    % Function to create a randome GSP network defined in an Ising
    % interaction matrix J and external field vector h. 

    J = zeros(numSpins); % Initialize random Spin matrix
    
    h = randn(numSpins,1);
    
    J(1,2) = randn();
%     J(1,3) = randn();
%     J(2,3) = randn();

    Jcons = (J ~= 0 ); % Use this to account for nodes we "dont" add. I.e phantom nodes

    for i = 3:numSpins
        [row, col] = find(Jcons);
        randidx = randi(length(row));

        if randi(2) == 2 % Connect two nodes
            J(row(randidx),i) = randn();
            J(col(randidx),i) = randn(); 
            Jcons(row(randidx),i) = 1;
            Jcons(col(randidx),i) = 1; 
        else % Connect to one
            if randi(2) == 2
                J(row(randidx),i) = randn();
                Jcons(row(randidx),i) = 1;
                Jcons(col(randidx),i) = 1; 
            else
                J(col(randidx),i) = randn(); 
                Jcons(row(randidx),i) = 1;
                Jcons(col(randidx),i) = 1; 
            end
        end
    end

    J = J + J.';

end

function [h1, J12, J13] = inverse_Ising_GSP_01_helper(m, C, step_size)
    % Inputs: Magnetizations m and correlations of 3-node system with 0,1
    % variables, where node 1 is the node we want to add back and nodes 2 and 3
    % are the neighbors of node 1. We also take as input the step size for the
    % nonlinear equation solver.
    %
    % Outputs: Effective external field on node 1 h1, and effective
    % interactions between node 1 and node 2 J12 and node 1 and node 3 J13.
    % These effective parameters are computed numerically using Eqs. 22-24 from
    % Report 1.
    %
    % NOTE: The difference between this function and 'inverse_Ising_GSP_helper'
    % is that here we consider an Ising system with 0,1 variables.
    
    % Threshold:
    threshold = 10^(-10);
    
    % Quantities of interest:
    m1 = m(1);
    m2 = m(2);
    m3 = m(3);
    C12 = C(1,2);
    C13 = C(1,3);
    C23 = C(2,3);
    
    % Initial guess at effective parameters (comes from the approximate
    % solution for small h1, J12, and J13):
    h1 = -2*((1 - 2*m1)*C23^2 - m2*m3*(1 + 2*C12 + 2*C13 - 2*m1 - m2 - m3) + 2*C23*(C13*m2 + (C12 - m2)*m3))...
        /(C23^2 - 2*C23*m2*m3 - m2*m3*(1 - m2 - m3));
    J12 = 4*(m1*m3*(m2 - C23) - C12*m3*(1 - m3) + C13*(C23 - m2*m3))...
        /(C23^2 - 2*m2*m3*C23 - m2*m3*(1 - m2 - m3));
    J13 = 4*(m1*m2*(m3 - C23) - C13*m2*(1 - m2) + C12*(C23 - m2*m3))...
        /(C23^2 - 2*m2*m3*C23 - m2*m3*(1 - m2 - m3));
    
    % Loop until convergence:
    diff = 1;
    count = 1;
    
    while diff > threshold
        
        % Compute magnetization and correlations for node i:
        m1_temp = (1 - m2 - m3 + C23)/(1 + exp(-h1)) + (m2 - C23)/(1 + exp(-J12 - h1))...
            + (m3 - C23)/(1 + exp(-J13 - h1)) + C23/(1 + exp(-J12 - J13 - h1));
        C12_temp = (m2 - C23)/(1 + exp(-J12 - h1)) + C23/(1 + exp(-J12 - J13 - h1));
        C13_temp = (m3 - C23)/(1 + exp(-J13 - h1)) + C23/(1 + exp(-J12 - J13 - h1));
        
        % Derivatives of magnetization and correlations with respect to parameters:
        dm1_dh1 = (1 - m2 - m3 + C23)*exp(-h1)/(1 + exp(-h1))^2 + (m2 - C23)*exp(-J12 - h1)/(1 + exp(-J12 - h1))^2 ...
            + (m3 - C23)*exp(-J13 - h1)/(1 + exp(-J13 - h1))^2 + C23*exp(-J12 - J13 - h1)/(1 + exp(-J12 - J13 - h1))^2;
        dm1_dJ12 = (m2 - C23)*exp(-J12 - h1)/(1 + exp(-J12 - h1))^2 ...
            + C23*exp(-J12 - J13 - h1)/(1 + exp(-J12 - J13 - h1))^2;
        dm1_dJ13 = (m3 - C23)*exp(-J13 - h1)/(1 + exp(-J13 - h1))^2 ...
            + C23*exp(-J12 - J13 - h1)/(1 + exp(-J12 - J13 - h1))^2;
        dC12_dh1 = dm1_dJ12;
        dC12_dJ12 = dm1_dJ12;
        dC12_dJ13 = C23*exp(-J12 - J13 - h1)/(1 + exp(-J12 - J13 - h1))^2;
        dC13_dh1 = dm1_dJ13;
        dC13_dJ12 = dC12_dJ13;
        dC13_dJ13 = dm1_dJ13;
        
        % Compute change in parameters:
        dhJ = [dm1_dh1, dm1_dJ12, dm1_dJ13; dC12_dh1, dC12_dJ12, dC12_dJ13;...
            dC13_dh1, dC13_dJ12, dC13_dJ13]\[m1_temp - m1; C12_temp - C12; C13_temp - C13];
        
        diff = sum(abs(dhJ));
        
        % Compute new parameters:
        h1 = h1 - step_size*dhJ(1);
        J12 = J12 - step_size*dhJ(2);
        J13 = J13 - step_size*dhJ(3); 
        
        % Check for nan values:
        if ((isnan(h1) + isnan(J12) + isnan(J13)) > 0) || max(abs([h1, J12, J13])) > 100
            return;
        end
        
        count = count + 1;
        if count > 1000
            disp('Did not finish')
            return;
        end
    end
end
    
function [entropy] = BoltzmanEnt(J, h, kT)
    %  USING 0,1 FOR SPINS

    numSpins = length(h);
   
    
    z = 0; % Partition Function
    probs = zeros(2^numSpins, 1);
    for i = 0:2^numSpins -1
    
        spin_set = decimalToBinaryVector(i, numSpins).';
        z = z + exp(-energyIsing(spin_set, J, h)/kT);

    end 
    
    entropy = 0;
    for i = 1:2^numSpins
        spin_set = decimalToBinaryVector(i -1, numSpins).';
        
        probs(i) =  exp(-energyIsing(spin_set, J, h)/kT)/z;

        entropy = entropy -  probs(i)*log(probs(i));
    end

    
    
end

function [meanspin, corr,threepoint] = Exact_Ising(J, h, kT)
    %  USING 0,1 FOR SPINS

    numSpins = length(h);
   
    
    z = 0; % Partition Function
    probs = zeros(2^numSpins, 1);
    for i = 0:2^numSpins -1
    
        spin_set = decimalToBinaryVector(i, numSpins).';
        z = z + exp(-energyIsing(spin_set, J, h)/kT);

    end 
    
    meanspin = zeros(numSpins, 1);
    corr = zeros(numSpins, numSpins);
    threepoint = zeros(numSpins, numSpins,numSpins);
    for i = 1:2^numSpins
        spin_set = decimalToBinaryVector(i -1, numSpins).';
        
        probs(i) =  exp(-energyIsing(spin_set, J, h)/kT)/z;
        meanspin = meanspin + spin_set.*probs(i);
        corr =  corr + probs(i)*(spin_set)*spin_set.';

        
        % Create a 3D grid of matrices from the input vectors
        [AA, BB, CC] = meshgrid(spin_set, spin_set, spin_set);
        
        % Multiply the corresponding elements element-wise
        tensor = AA .* BB .* CC;


        threepoint = threepoint + probs(i)*tensor;
    end 
    
    
end

function Emean = energyIsing(spin, J, h)
    %ENERGYISING Mean energy per spin.
    %   Emean = ENERGYISING(spin, J) returns the mean energy per spin of the
    %   configuration |spin|. |spin| is a matrix of +/- 1's. |J| is a scalar.

    Emean = -(1./2).*(spin.' * J * spin) + trace(J)  - h.'*spin;
end
