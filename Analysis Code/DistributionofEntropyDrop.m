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
%%

hind = log(datamean);

datacorr_pseudo = datacorr_pseudo - eye(num_nuerons)/(num_bins +1);

Jpair = log(datacorr_pseudo) + log(1 + datacorr_pseudo - datamean - datamean.') - log(datamean - datacorr_pseudo) - log(datamean.' - datacorr_pseudo);

Jpair = real(Jpair);

%%
numSpins = length(datamean);

spins = 1:numSpins;

uniquePairs = zeros(numSpins*(numSpins-1)/2, 2);
for i = 1:numSpins-1
    for j = i+1:numSpins
        uniquePairs( (i - 1)*numSpins - (i -1)*i/2 + j-i, 1) = i;
        uniquePairs((i - 1)*numSpins - (i -1)*i/2 + j-i, 2) = j;
    end
end


%%

[~, randPerms] = sort(rand(size(uniquePairs, 1), 1));

%%

deltaHs = zeros(numSpins, 100);

tic
for i = 1:10000
    

    
    deltaH = pointthreeMutual(Jpair, hind, datamean, datacorr_pseudo, numSpins, 1:numSpins, uniquePairs(randPerms(i),:));
    
    
    deltaHs(:,i) = deltaH.';

end
toc
%%

deltaHs = sort(reshape(deltaHs,[],1), 'descend');

loglog(1:length(deltaHs),deltaHs, '.')
xlabel('Rank')
ylabel('Entropy Drop GSP')

function [deltaH]  = pointthreeMutual(J, h, datamean, datacorr, numSpins, spins, pair)
    % Calculate the DKL drop for a set of spins connecting to a pair of
    % spins

    deltaH = zeros(numSpins,1);
    j = pair(1);
    k = pair(2);
    
    h(j) = log(datamean(j) - datacorr(j,k)) - log(1+datacorr(j,k) - datamean(j) - datamean(k));
    h(k) = log(datamean(k) - datacorr(j,k)) - log(1+datacorr(j,k) - datamean(j) - datamean(k));
    
    
    for i = spins
        
        if i == j || i == k
            deltaH(i) = 0; 
            continue
        end

        % ----- %
        m = datamean([i,j,k]);
        C = datacorr([i,j,k],[i,j,k]);
        
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


        % ---- %

        deltaH(i) = -log(exp(h(i))+1) + Jij*datacorr(i,j) + Jik*datacorr(i,k) + h(i)*datamean(i) - datamean(i)*log(datamean(i))- (1-datamean(i))*log(1-datamean(i)) + (Jjk - Jjkold)*datacorr(j,k) + (hj - hjold)*datamean(j) + (hk- hkold)*datamean(k); 
    end

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
    

function deltaS = MI2(mean,corr)
    
    a = -mean.*log(mean) - (1-mean).*log(1-mean) - (mean.*log(mean)).' - ((1-mean).*log(1-mean)).';
    b = log(corr.^corr); % Hack to make 0*log 0 = 0
    c = (mean - corr).*log(mean - corr) + ((mean - corr).*log(mean - corr)).';
    d = (1 - mean - mean.' +corr).*log((1 - mean - mean.' +corr));
    deltaS = real(a + b + c + d);

    for i = 1:length(mean)
        deltaS(i,i) = 0; %-mean(i)*log(mean(i)) - (1-mean(i))*log(1-mean(i));
    end

    deltaS(isnan(deltaS)) = 0;
end
