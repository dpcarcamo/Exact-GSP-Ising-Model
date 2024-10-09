%% Frustration Analysis
clear
clc

currentPath = pwd;
folders = split(currentPath, '\');
newPath = join(folders(1:length(folders)-2),"\");

addpath(strcat(newPath{1} , '\Model Data'))

listing = dir(strcat(newPath{1} , '\Model Data'));

filename = listing(8).name


% Add Helper Function to Path
newPath = join(folders(1:length(folders)-1),"\");
addpath(strcat(newPath{1}, '\Helper Function'))
load(filename)

%% Optimal

[G,D] = decimate_GSP(JGSP);

frustrations = zeros(num_nuerons-2,1);

count = 1;
for keep = D(1:size(D,1)-1,:).'
    
    
    frustrations(count) = JGSP(keep(1),keep(2))*JGSP(keep(2),keep(3))*JGSP(keep(3),keep(1));

    count = count + 1;

    
end

histogram(frustrations, 'Normalization','probability',EdgeColor='none', FaceAlpha=0.7)
hold on

% Random

[G,D] = decimate_GSP(JGSPrand);

frustrations = zeros(num_nuerons-2,1);

count = 1;
for keep = D(1:size(D,1)-1,:).'
    
    
    frustrations(count) = JGSPrand(keep(1),keep(2))*JGSPrand(keep(2),keep(3))*JGSPrand(keep(3),keep(1));

    count = count + 1;

    
end

histogram(frustrations, 'Normalization','probability', EdgeColor='none', FaceAlpha=0.8)

% Min Dist

[G,D] = decimate_GSP(JGSPdist);

frustrations = zeros(num_nuerons-2,1);

count = 1;
for keep = D(1:size(D,1)-1,:).'
    
    
    frustrations(count) = JGSPdist(keep(1),keep(2))*JGSPdist(keep(2),keep(3))*JGSPdist(keep(3),keep(1));

    count = count + 1;

    
end

histogram(frustrations, 'Normalization','probability', EdgeColor='none', FaceAlpha=0.6)

hold off
legend('Optimal','Random','Min Dist')
xlim([-15,15])
xlabel('JijJikJjk')
ylabel('Probability')


%% Synergy/Redunecy 

J = JGSP;
[G,D] = decimate_GSP(J);

frustrations = zeros(num_nuerons-2,1);

Hind = Entropy(datamean);

Mutual_info = MI2(datamean,datacorr_pseudo);

synergy = zeros(num_nuerons-2,1);
count = 1;
for keep = D(1:size(D,1)-1,:).'
    
    
    frustrations(count) = J(keep(1),keep(2))*J(keep(2),keep(3))*J(keep(3),keep(1));

    threepoint = sum(binary(:,keep(1)).*(binary(:,keep(2)).*binary(:,keep(3))))/6956;

    Hthree = Entropy_Three(datamean(keep), datacorr_pseudo(keep,keep), threepoint);

    Hdrop = Hind(keep(1)) + Hind(keep(2)) + Hind(keep(3)) - Hthree - Mutual_info(keep(1),keep(2))  - Mutual_info(keep(1),keep(3))  - Mutual_info(keep(3),keep(2));


    synergy(count) = Hdrop;

    count = count + 1;
    
end

histogram(synergy, 'Normalization','probability', EdgeColor='none', FaceAlpha=0.6)
hold on


J = JGSPrand;
[G,D] = decimate_GSP(J);

frustrations = zeros(num_nuerons-2,1);

Hind = Entropy(datamean);

Mutual_info = MI2(datamean,datacorr_pseudo);

synergy = zeros(num_nuerons-2,1);
count = 1;
for keep = D(1:size(D,1)-1,:).'
    
    
    frustrations(count) = J(keep(1),keep(2))*J(keep(2),keep(3))*J(keep(3),keep(1));

    threepoint = sum(binary(:,keep(1)).*(binary(:,keep(2)).*binary(:,keep(3))))/6956;

    Hthree = Entropy_Three(datamean(keep), datacorr_pseudo(keep,keep), threepoint);

    Hdrop = Hind(keep(1)) + Hind(keep(2)) + Hind(keep(3)) - Hthree - Mutual_info(keep(1),keep(2))  - Mutual_info(keep(1),keep(3))  - Mutual_info(keep(3),keep(2));


    synergy(count) = Hdrop;

    count = count + 1;
    
end

histogram(synergy, 'Normalization','probability', EdgeColor='none', FaceAlpha=0.6)

J = JGSPdist;
[G,D] = decimate_GSP(J);

frustrations = zeros(num_nuerons-2,1);

Hind = Entropy(datamean);

Mutual_info = MI2(datamean,datacorr_pseudo);

synergy = zeros(num_nuerons-2,1);
count = 1;
for keep = D(1:size(D,1)-1,:).'
    
    
    frustrations(count) = J(keep(1),keep(2))*J(keep(2),keep(3))*J(keep(3),keep(1));

    threepoint = sum(binary(:,keep(1)).*(binary(:,keep(2)).*binary(:,keep(3))))/6956;

    Hthree = Entropy_Three(datamean(keep), datacorr_pseudo(keep,keep), threepoint);

    Hdrop = Hind(keep(1)) + Hind(keep(2)) + Hind(keep(3)) - Hthree - Mutual_info(keep(1),keep(2))  - Mutual_info(keep(1),keep(3))  - Mutual_info(keep(3),keep(2));


    synergy(count) = Hdrop;

    count = count + 1;
    
end

histogram(synergy, 'Normalization','probability', EdgeColor='none', FaceAlpha=0.6)




hold off
legend('Optimal','Random','Min Dist')

ylabel('Probability')
xlim([-0.03,0.01])

%%
figure
plot(frustrations, synergy, '.', MarkerSize=8, Color=[0 0.4470 0.7410])
hold on
plot(linspace(-100,75),zeros(100,1),  '--', 'Color',[0.4660 0.6740 0.1880], 'LineWidth',2)
plot(zeros(100,1),linspace(-0.08,0.02), '--', 'Color',[0.4660 0.6740 0.1880], 'LineWidth',2)
hold off
ylabel('Synergy')
xlabel('Frustration')

%%

left = synergy < 0;
up = frustrations > 0;

sum(left.*up)
sum((1-left).*(up))
sum((1-left).*(1-up))
sum((left).*(1-up))


%% Synergy/Redunecy using model fits

J = JGSP;
[G,D] = decimate_GSP(J);

frustrations = zeros(num_nuerons-2,1);

Hind = Entropy(datamean);

Mutual_info = MI2(datamean,datacorr_pseudo);

synergy = zeros(num_nuerons-2,1);
count = 1;
for keep = D(1:size(D,1)-1,:).'
    
    [h, J, Hthree, ~] = find_GSP_update(datamean(keep),datacorr_pseudo(keep,keep));

    frustrations(count) = J(1,2)*J(2,3)*J(1,3);

    %threepoint = sum(binary(:,keep(1)).*(binary(:,keep(2)).*binary(:,keep(3))))/6956;

    %Hthreett = Entropy_Three(datamean(keep), datacorr_pseudo(keep,keep), threepoint);

    Hdrop = Hind(keep(1)) + Hind(keep(2)) + Hind(keep(3)) - Hthree - Mutual_info(keep(1),keep(2))  - Mutual_info(keep(1),keep(3))  - Mutual_info(keep(3),keep(2));


    synergy(count) = Hdrop;

    count = count + 1;
    
end

% histogram(synergy, 'Normalization','probability', EdgeColor='none', FaceAlpha=0.6)
% 
% hold off
% legend('Optimal')
% 
% ylabel('Probability')

plot(frustrations, synergy, '.', MarkerSize=8, Color=[0 0.4470 0.7410])
hold on
plot(linspace(-150,100),zeros(100,1),  '--', 'Color',[0.4660 0.6740 0.1880], 'LineWidth',2)
plot(zeros(100,1),linspace(-0.08,0.02), '--', 'Color',[0.4660 0.6740 0.1880], 'LineWidth',2)

ylabel('Synergy')
xlabel('Frustration')

%%

left = synergy < 0;
up = frustrations > 0;

sum(left.*up)
sum((1-left).*(up))
sum((1-left).*(1-up))
sum((left).*(1-up))

%% Comparision Between J true and J model

J = JGSP;
[G,D] = decimate_GSP(J);

Jstrue = zeros(3,num_nuerons-2);
Jsmodel = zeros(3,num_nuerons-2);

frustrations = zeros(num_nuerons-2,1);

Hind = Entropy(datamean);

Mutual_info = MI2(datamean,datacorr_pseudo);

synergy = zeros(num_nuerons-2,1);
count = 1;
for keep = D(1:size(D,1)-1,:).'
    
    [h, Jtemp, Hthree, ~] = find_GSP_update(datamean(keep),datacorr_pseudo(keep,keep));

    Jstrue(:,count) = [Jtemp(1,2),Jtemp(1,3),Jtemp(2,3)];
    Jsmodel(:,count) = [J(keep(1),keep(2)),J(keep(1),keep(3)),J(keep(2),keep(3))];


    count = count + 1;
    
end


%%
figure
plot(Jstrue.',Jsmodel.', '.', Color=[0, 0.4470,0.7410])
xlabel('J values from Triangle Fit')
ylabel('J values from Full Fit')

%% Synergy/Redunecy sweep


frustrations = zeros(10^7,1);


synergy = zeros(10^7,1);


for i = 1:10^6

    Jxyxyz = 5*randn(1);
    Jxzxyz = 5*randn(1);
    Jyzxyz = 5*randn(1);

    hxxyz = 2*randn(1);
    hyxyz = 2*randn(1);
    hzxyz = 2*randn(1);

    frustrations(i) = Jxyxyz*Jxzxyz*Jyzxyz;


    Jxyxy = Jxyxyz + log(exp(hzxyz) + 1) - log(exp(hzxyz + Jxzxyz) + 1) - log(exp(hzxyz + Jyzxyz) + 1) + log(exp(hzxyz + Jxzxyz + Jyzxyz) + 1);
    Jxzxz = Jxzxyz + log(exp(hyxyz) + 1) - log(exp(hyxyz + Jyzxyz) + 1) - log(exp(hyxyz + Jxyxyz) + 1) + log(exp(hyxyz + Jyzxyz + Jxyxyz) + 1);
    Jyzyz = Jyzxyz + log(exp(hxxyz) + 1) - log(exp(hxxyz + Jxzxyz) + 1) - log(exp(hxxyz + Jxyxyz) + 1) + log(exp(hxxyz + Jxzxyz + Jxyxyz) + 1);
    
    hxxy = hxxyz - log(exp(hzxyz) + 1) + log(exp(hzxyz + Jxzxyz) + 1);
    hxxz = hxxyz - log(exp(hyxyz) + 1) + log(exp(hyxyz + Jxyxyz) + 1);
    hyxy = hyxyz - log(exp(hzxyz) + 1) + log(exp(hzxyz + Jyzxyz) + 1);
    hyyz = hyxyz - log(exp(hxxyz) + 1) + log(exp(hxxyz + Jxyxyz) + 1);
    hzxz = hzxyz - log(exp(hyxyz) + 1) + log(exp(hyxyz + Jyzxyz) + 1);
    hzyz = hzxyz - log(exp(hxxyz) + 1) + log(exp(hxxyz + Jxzxyz) + 1);
    
    hxx = hxxy - log(exp(hyxy) + 1) + log(exp(hyxy + Jxyxy) + 1);
    hyy = hyxy - log(exp(hxxy) + 1) + log(exp(hxxy + Jxyxy) + 1);
    hzz = hzxz - log(exp(hxxz) + 1) + log(exp(hxxz + Jxzxz) + 1);
    
    fxy = log(exp(hzxyz) + 1);
    fxz = log(exp(hyxyz) + 1);
    fyz = log(exp(hxxyz) + 1);
    
    fx = fxy + log(exp(hyxy) + 1);
    fy = fxy + log(exp(hxxy) + 1);
    fz = fxz + log(exp(hxxz) + 1);
    
    f = fz + log(exp(hzz) + 1);
    
    Zxyz = exp(f) ;
    Zxy = Zxyz*exp(-fxy);
    Zxz = Zxyz*exp(-fxz);
    Zyz = Zxyz*exp(-fyz);
    
    Zx = Zxyz*exp(-fx);
    Zy = Zxyz*exp(-fy);
    Zz = Zxyz*exp(-fz);

    Syn = log((Zxy *Zxz* Zyz)/(Zxyz *Zx* Zy* Zz)) - exp(hxx)* (hxxy + hxxz - hxxyz - hxx)/Zx - exp(hyy) *(hyxy + hyyz - hyxyz - hyy)/Zy - exp(hzz) *(hzyz + hzxz - hzxyz - hzz)/Zz - exp(hxxy + hyxy + Jxyxy) *(Jxyxy - Jxyxyz)/Zxy - exp(hxxz + hzxz + Jxzxz) *(Jxzxz - Jxzxyz)/Zxz  - exp(hyyz + hzyz + Jyzyz)* (Jyzyz - Jyzxyz)/Zyz;

    synergy(i) = Syn;

end

%%
figure
plot(frustrations, synergy, '.', MarkerSize=8, Color=[0 0.4470 0.7410])
hold on
plot(linspace(-50,50),zeros(100,1),  '--', 'Color',[0.4660 0.6740 0.1880], 'LineWidth',2)
plot(zeros(100,1),linspace(-0.3,0.2), '--', 'Color',[0.4660 0.6740 0.1880], 'LineWidth',2)
hold off
ylabel('Synergy')
xlabel('Frustration')

left = synergy < 0;
up = frustrations > 0;

sum(left.*up)
sum((1-left).*(up))
sum((1-left).*(1-up))
sum((left).*(1-up))

%%

frusts = linspace(-30,20,200);
boundry = zeros(size(frusts));
count = 1;

for i = frusts
    
    %bin = frustrations > i;

    boundry(count) = max(synergy(frustrations > i));
    count = count + 1;
end


figure
plot(frusts, boundry, '.')
xlabel('log(Frustration)')
ylabel('log(Synergy Boundry)')


%%
clc

data = csvread('data.csv');
%%
datatemp = flip(data);

boundry = zeros(size(data));
prev = 0;
for iter = 1:size(data,1)

    if prev - datatemp(iter,2) < 10^(-4)
        prev = datatemp(iter,2);
        boundry(iter,:) = datatemp(iter,:);

    end

end

boundry(boundry(:,2)==0,:) = [];


%plot(boundry(:,1),boundry(:,2), '.')

xShade = boundry(:,1);
yShade = boundry(:,2);
%plot(frustrations, synergy, '.', MarkerSize=8, Color=[0 0.4470 0.7410])
plotvariance(frustrations, synergy,15,'b')
hold on
%plot(linspace(-150,100),zeros(100,1),  '--', 'Color',[0.4660 0.6740 0.1880], 'LineWidth',2)
%%plot(zeros(100,1),linspace(-0.08,0.02), '--', 'Color',[0.4660 0.6740 0.1880], 'LineWidth',2)

ylabel('Synergy')
xlabel('Frustration')
fill([xShade.',-100,100], [yShade.',1,1], 'r', 'FaceAlpha', 0.3, 'EdgeAlpha',0);
hold off
ylim([-0.1,0.2])
xlim([min(xShade),max(xShade)])


function [H] = Entropy_Three(mean, corr, threepoint)


    xi = 1;
    xj = 2;
    xk = 3;
    Probs = zeros(2,2,2);

    Probs(2,2,2) = threepoint;
    Probs(1,2,2) = corr(xj,xk) - Probs(2,2,2);
    Probs(2,1,2) = corr(xi,xk) - Probs(2,2,2);
    Probs(1,1,2) = mean(xk) - corr(xi,xk) - corr(xj,xk) + Probs(2,2,2);

    Probs(2,2,1) = corr(xi,xj) - Probs(2,2,2);
    Probs(1,2,1) = mean(xj) - corr(xi,xj) - corr(xj,xk) + Probs(2,2,2);
    Probs(2,1,1) = mean(xi) - corr(xi,xj) - corr(xi,xk) + Probs(2,2,2);
    Probs(1,1,1) = 1 - Probs(2,2,2) - Probs(2,2,1) - Probs(2,1,2) - Probs(1,2,2) - Probs(2,1,1) - Probs(1,2,1) - Probs(1,1,2) ;


    H = real(sum(-Probs.*log(Probs), 'all', "omitnan"));

end
