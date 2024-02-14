%% Frustration Analysis
clear
clc

currentPath = pwd;
folders = split(currentPath, '\');
newPath = join(folders(1:length(folders)-2),"\");

addpath(strcat(newPath{1} , '\Model Data'))

listing = dir(strcat(newPath{1} , '\Model Data'));

filename = listing(8).name

load(filename)
% Add Helper Function to Path
newPath = join(folders(1:length(folders)-1),"\");
addpath(strcat(newPath{1}, '\Helper Function'))


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

% histogram(synergy, 'Normalization','probability', EdgeColor='none', FaceAlpha=0.6)
% 
% hold off
% legend('Optimal')
% 
% ylabel('Probability')

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


for i = 1:10^7

    Jxyxyz = randn(1);
    Jxzxyz = randn(1);
    Jyzxyz = randn(1);

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

%
figure
plot(frustrations, synergy, '.', MarkerSize=8, Color=[0 0.4470 0.7410])
hold on
plot(linspace(-200,200),zeros(100,1),  '--', 'Color',[0.4660 0.6740 0.1880], 'LineWidth',2)
plot(zeros(100,1),linspace(-0.6,0.2), '--', 'Color',[0.4660 0.6740 0.1880], 'LineWidth',2)
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

frusts = -logspace(-6,1.2,100);
boundry = zeros(size(frusts));
count = 1;

for i = frusts
    
    %bin = frustrations > i;

    boundry(count) = max(synergy(frustrations > i));
    count = count + 1;
end


figure
plot(log10(-frusts), log10(boundry), '.')
xlabel('log(Frustration)')
ylabel('log(Synergy Boundry)')



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
