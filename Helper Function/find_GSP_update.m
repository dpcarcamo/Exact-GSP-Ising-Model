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
        count = count + 1;
    end

end
