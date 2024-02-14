% Script to analyze all different iterations of networks.
clear
clc

currentPath = pwd;
folders = split(currentPath, '\');
newPath = join(folders(1:length(folders)-2),"\");

listing = dir(strcat(newPath{1} , '\Stringer Data\Natimg2800'));
addpath(strcat(newPath{1} , '\Stringer Data\Natimg2800'))

% Add Helper Function to Path
newPath = join(folders(1:length(folders)-1),"\");
addpath(strcat(newPath{1}, '\Helper Function'))


for namenum = 3:9

    filename = listing(namenum).name;
    
    ExpData = matfile(filename);
    
    
    %
    Stim = ExpData.stim;
    Response = Stim.resp;
    Spont = Stim.spont;
    Locations = ExpData.med;
    dbinfo = ExpData.db;
    mouse_name = dbinfo.mouse_name;
    
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
    
    
    
    
    % Distance
    
    distanceSMetric = squareform(pdist(Locations));
    distanceMetric = distanceSMetric(5006,:);
    
    %
    
    
    [sortDistances, I] = sort(distanceMetric);
    
    
    %
    randtrial = 50;
    Data = zeros(20,num_nuerons, 45);  
      
    
    count = 1;
    for j = unique([round(logspace(log10(2),log10(num_nuerons/500), 10))  , round(logspace(log10(num_nuerons/500),log10(num_nuerons/100), 4))])
        i = count;

        centerspins = [];
        for k = 1:round(num_nuerons/j)
            
            spins = 1:num_nuerons;
            spins(centerspins) = [];
            randspin = spins(randi(length(spins)));
            centerspins = [centerspins, randspin];
            

            distanceMetric = distanceSMetric(randspin,:);
            [sortDistances, I] = sort(distanceMetric);
            region   = (distanceMetric <=  sortDistances(j));
            
            reduceddatamean = datamean(region);
            reduceddatacorr = datacorr_pseudo(region, region);
            reducedLocation = Locations(region,:);
            
            Hind = sum(Entropy(reduceddatamean));
            
            
            % Minimal distance tree
    
            [JTreedist, hTreedist, HTreedist] = Minimum_Distance_Tree(reduceddatamean,reduceddatacorr,reducedLocation);
    
            % Minimal distance GSP
    
            [JGSPdist, hGSPdist, HGSPdist] = Minimum_Distance_GSP(reduceddatamean,reduceddatacorr,reducedLocation);
    
            % Find Optimal Tree
    
            [JTree, hTree, HTree] = FindTree(reduceddatamean, reduceddatacorr);
            
            avHGSPrand = 0;
            avmeanhGSPrand = 0;
            avmeanhintGSPrand = 0;
            avTreeRand = 0;
            avmeanhTreeRand = 0;
            avmeanhintTreeRand = 0;
    
            for ra = 1:randtrial
                % Random GSP
                
                [JGSPrand, hGSPrand, HGSPrand] = RandomGSPFit(reduceddatamean, reduceddatacorr);
        
                avHGSPrand = avHGSPrand + HGSPrand;
                avmeanhGSPrand = avmeanhGSPrand + mean(hGSPrand/2 + sum(JGSPrand));
                avmeanhintGSPrand = avmeanhintGSPrand + mean(JGSPrand*(2*reduceddatamean-1).'/4);
        
                
                %Random Tree
                
                [JTreeRand, hTreeRand, HTreeRand] = Random_Tree(reduceddatamean, reduceddatacorr);
    
                avTreeRand = avTreeRand + HTreeRand;
                avmeanhTreeRand = avmeanhTreeRand + mean(hTreeRand/2 + sum(JTreeRand));
                avmeanhintTreeRand = avmeanhintTreeRand + mean(JTreeRand*(2*reduceddatamean-1).'/4);
            end
        
            avHGSPrand = avHGSPrand/randtrial;
            avmeanhGSPrand = avmeanhGSPrand/randtrial;
            avmeanhintGSPrand = avmeanhintGSPrand/randtrial;
            avTreeRand = avTreeRand/randtrial;
            avmeanhTreeRand = avmeanhTreeRand/randtrial;
            avmeanhintTreeRand = avmeanhintTreeRand/randtrial;
            
            % Find Optimal GSP
            
            
            [hGSP, JGSP, HGSP] = find_GSP_update(reduceddatamean, reduceddatacorr);
            
    
            %Entropy Drops
            
            Data(1, k, i) = Hind;
            Data(2, k, i) = HGSPdist;
            Data(3, k, i) = HTreedist;
            Data(4, k, i) = avHGSPrand;
            Data(5, k, i) = avTreeRand;
            Data(6, k, i) = HGSP;
            Data(7, k, i) = HTree;
            
            Data(8, k, i) = mean(hGSPdist/2+sum(JGSPdist));
            Data(9, k, i) = mean(hTreedist/2 + sum(JTreedist));
            Data(10, k, i) = avmeanhGSPrand;
            Data(11, k, i) = avmeanhTreeRand;
            Data(12, k, i) = mean(hGSP/2 + sum(JGSP));
            Data(13, k, i) = mean(hTree/2 + sum(JTree));
    
            Data(14, k, i) = mean(JGSPdist*(2*reduceddatamean-1).'/4);
            Data(15, k, i) = mean(JTreedist*(2*reduceddatamean-1).'/4);
            Data(16, k, i) = avmeanhintGSPrand;
            Data(17, k, i) = avmeanhintTreeRand;
            Data(18, k, i) = mean(JGSP*(2*reduceddatamean-1).'/4);
            Data(19, k, i) = mean(JTree*(2*reduceddatamean-1).'/4);
            Data(20, k, i) = j;
        end
        count = count + 1
        j
    end
    
    save(strcat(filename(1:end-4), '_fitted'))
end