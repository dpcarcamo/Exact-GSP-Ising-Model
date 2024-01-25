clear all
clc
numSpins = 10;

trials = 10;
Data = zeros(trials,10);

tic
for i = 1:trials


    [J,h] = RandomGSPIsing(numSpins);
    
    %[ExactMean, ExactCorr, threepoint] = Exact_Ising(J, h, 1); % Brute Force
    [exactGSPmeanT, exactGSPcorrT] = correlations_GSP_01(J,h);
    H = Entropy_GSP(J,h, exactGSPmeanT, exactGSPcorrT);
    
    %[hguess, Jguess, EntModelGSP, NewEnt] = find_GSP(exactGSPmeanT, exactGSPcorrT, threepoint);
    [hguess2, Jguess2, EntModelGSP2, NewEnt2] = find_GSP_update(exactGSPmeanT, exactGSPcorrT);

    [exactGSPmeanG, exactGSPcorrG] = correlations_GSP_01(Jguess2,hguess2);
    
    

    bJ = (J ~= 0);
    %bJguess = (round(Jguess, 7) ~= 0 );
    bJguess2 = (round(Jguess2, 7) ~= 0 );

    [JT, HT, EntModelTree] = FindTree(exactGSPmeanT, exactGSPcorrT);

    bJT = (JT ~= 0);
    sum(bJ(bJT == 1), "all")/2;
% 
%     Data(i,1) = sum(abs(h-hguess), "all")/numSpins;
%     Data(i,2) = sum(abs(J-Jguess), "all")/(numSpins.*(numSpins-1)/2);
%     diffcons = bJ - bJguess;
%     Data(i,3) = sum(diffcons(diffcons > 0), "all")/2; % Missing connections
%     Data(i,4) = -sum(diffcons(diffcons < 0), "all")/2; % Over connected
%     Data(i,5) = sum(abs(exactGSPcorrG-exactGSPcorrT -exactGSPmeanG+exactGSPmeanT)/(numSpins.*(numSpins-1)),"all");
%     Data(i,6) = sum(abs(exactGSPmeanG-exactGSPmeanT)/numSpins,"all");
%     Data(i,7) =  (numSpins-1 - sum(bJ(bJT == 1), "all")/2)/(numSpins-1); % Compare to best tree
%     Data(i,8) = EntModelTree;
%     Data(i,9) = NewEnt;
%     Data(i,10) = H;
%     Data(i, 11) = sum(ExactMean - exactGSPmeanG, 'all');
%     Data(i, 12) = sum(ExactCorr -exactGSPcorrG, "all");
%     Data(i,13) = EntModelGSP;

    %diffcons = bJ - bJguess;
    %Data(i,1) = sum(diffcons(diffcons > 0), "all")/2; % Missing connections
    %Data(i,2) = -sum(diffcons(diffcons < 0), "all")/2; % Over connected
    diffcons2 = bJ - bJguess2;
    Data(i,3) = sum(diffcons2(diffcons2 > 0), "all")/2; % Missing connections
    Data(i,4) = -sum(diffcons2(diffcons2 < 0), "all")/2; % Over connected
    Data(i,5) = H;
    %Data(i,6) = BoltzmanEnt(J, h, 1);
    %Data(i,7) = EntModelGSP;
    %Data(i,8) = NewEnt;
    Data(i,9) = EntModelGSP2;
    Data(i,10) = NewEnt2;
    Data(i, 11) = Entropy_GSP(Jguess2, hguess2, exactGSPmeanG, exactGSPcorrG);
    %Data(i,11) = BoltzmanEnt(temp1, hguess, 1);
    %Data(i,12) = BoltzmanEnt(temp, hguess2, 1);


end
toc
disp('Done')
disp(mean(Data))

%% Timing 
clear
clc
trials = [4,10,20,50,100,150,200,500, 1000];

times = zeros(length(trials), 1);

count = 1;
for numSpins = trials       
    
    [J,h] = RandomGSPIsing(numSpins);
    disp('Generated Network')
    [exactGSPmeanT, exactGSPcorrT] = correlations_GSP_01(J,h);
    disp('Found Stats')
    tStart = tic;
    [hmodel, Jmodel, Entmodel, NewEnt] = find_GSP_update(exactGSPmeanT, exactGSPcorrT);
    times(count) = toc(tStart);
    count = count + 1
end


%%
loglog(trials,times, 'o')
hold on
xlabel('Size (number of nodes)', 'FontSize',14)
ylabel('Time to find network (seconds)', 'FontSize',14)
c = polyfit(log10(trials(3:end)),log10(times(3:end)),1)
loglog(trials,(10^c(2))*trials.^c(1))
legend('Data', append('Slope = ' , num2str(c(1))), 'Location','northwest', 'Fontsize',14)
hold off
%%
clc
numSpins = 2;
tic
[J,h] = RandomGSPIsing(numSpins);
toc

tic
[ExactMean, ExactCorr, threepoint] = Exact_Ising(J, h, 1);
toc
%%
H = Entropy(ExactMean);

G = graph(J, 'upper');


% Name nodes with their indices
nodeNames = cellstr(num2str((1:numSpins).')); % Convert indices to cell array of strings
G.Nodes.Name = nodeNames; % Assign node names


plot(G)
hold on
tic
[hguess, Jguess, EntModelGSP, NewEnt] = find_GSP_update(ExactMean, ExactCorr);
toc
tic
[hguess2, Jguess2, EntModelGSP2, NewEnt2] = find_GSP_update(ExactMean, ExactCorr);
toc
plot(graph(round(Jguess2,7)))
[JT, HT, EntModelTree] = FindTree(ExactMean, ExactCorr);

EntModelGSP, NewEnt
EntModelGSP2, NewEnt2, EntModelTree, sum(H)

hguess- h;
Jguess- J;
EntModelGSP - H;



%plplot(Gguess)
%nodeNames = cellstr(num2str((1:numSpins).')).';
%Jguess = full(adjacency(reordernodes(Gguess,nodeNames), 'weighted'));
H = graph(round(Jguess, 7));
%plot(H)


% [Jguess, hguess, EntModelTree] = FindTree(ExactMean, ExactCorr);
% EntModelTree
H = Entropy_GSP(J,h, ExactMean, ExactCorr);

% disp(EntModelTree - H)
% 
% plot(graph(Jguess))
hold off

% tic
% [PredictMean, PredictCorr, Predictthreepoint] = Exact_Ising(Jguess, hguess, 1);
% toc
% 
% [m, C, X, Z] = correlations_GSP_01(Jguess, hguess);
% 
% dmean = PredictMean - ExactMean;
% dcorr =PredictCorr -  ExactCorr;
% dthree = Predictthreepoint - threepoint;
% dcorr(Jguess > 0);
% [row, col] = find(D == 0);
% [G, D] = decimate_GSP(Jguess)
% D(row,:) = ones(size(D(row,:)));
% dthree(D);
%%


[G,D] = decimate_GSP(J)
Entropy_GSP(J,h, ExactMean, ExactCorr)
BoltzmanEnt(J, h, 1)

%%
clc
numSpins = 3;
tic
[J,h] = RandomGSPIsing(numSpins);
toc


% T = 1;
% f = T*log(exp(h/T)+1);
% df = 0 + log(exp(h/T)+1) + 1/(exp(-h/T)+1)*(0 - h/T);
% ddf = 0 + 1/(exp(-h/T)+1)*0 + exp(-h/T)/(T*(exp(-h/T)+1)^2)*((0 - h/T)^2);
% 
% -h/(1+exp(-h/T))
% -f + T*df
% 
% (h^2)*exp(-h/T)/(1+exp(-h/T))^2
% T*ddf


%
tic
[ExactMean, ExactCorr, threepoint, C] = Exact_Ising(J, h, 2);
toc

Emean = -(1./2).*sum(J.*ExactCorr,'all')  - h.'*ExactMean;
[Energy, Specific] = Thermodynamics(J,h, 2);

Emean - Energy

C- Specific
%%

clc
clear
numSpins = 100;

[J,h] = RandomGSPIsing(numSpins);

%%
Ts = linspace(0.1,3);

Specifics = zeros(size(Ts));

count = 1;
for T = Ts
    

    % [ExactMean, ExactCorr, threepoint, C] = Exact_Ising(J, h, T);
    % [ExactMean, ExactCorr]  = correlations_GSP_01(J/T,h/T);
    
    
    % Emean = -(1./2).*sum(J.*ExactCorr,'all')  - h.'*ExactMean;
    [Energy, Specific] = Thermodynamics(JGSP,hGSP, T);
    Specifics(count) = Specific;
    count = count + 1;

end
plot(Ts,Specifics, 'LineWidth',2)
xlabel('Temperature (k_bT)', 'FontSize',18)
ylabel('Specific Heat', 'FontSize',18)