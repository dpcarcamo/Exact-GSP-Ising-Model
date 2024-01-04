clear
clc

currentPath = pwd;
folders = split(currentPath, '\');
newPath = join(folders(1:length(folders)-2),"\");

addpath(strcat(newPath{1} , '\Model Data'))

load("natimg2800_2017-08-20_Fitting.mat")


%%

metrics = zeros(num_nuerons,num_nuerons,10);
paths  = zeros(num_nuerons,num_nuerons,10);

bJ = JGSPrand ~= 0;
distanceMetric = squareform(pdist(Locations));
distanceMetric = distanceMetric.*bJ;
distanceMetric = double(distanceMetric);
distanceMetric(distanceMetric ==0) = Inf;
distanceMetric(1:1+size(distanceMetric,1):end) = 0;
[metric,path] = FloydWarshall(distanceMetric);
metrics(:,:,1) = metric;
paths(:,:,1) = path;

bJ = double(bJ);
bJ(bJ ==0) = Inf;
bJ(1:1+size(bJ,1):end) = 0;
[metric,path] = FloydWarshall(bJ);
metrics(:,:,2) = metric;
paths(:,:,2) = path;

bJ = JTreedist ~= 0;
distanceMetric = squareform(pdist(Locations));
distanceMetric = distanceMetric.*bJ;
distanceMetric = double(distanceMetric);
distanceMetric(distanceMetric ==0) = Inf;
distanceMetric(1:1+size(distanceMetric,1):end) = 0;
[metric,path] = FloydWarshall(distanceMetric);
metrics(:,:,3) = metric;
paths(:,:,3) = path;

bJ = double(bJ);
bJ(bJ ==0) = Inf;
bJ(1:1+size(bJ,1):end) = 0;
[metric,path] = FloydWarshall(bJ);
metrics(:,:,4) = metric;
paths(:,:,4) = path;

bJ = JGSPdist ~= 0;
distanceMetric = squareform(pdist(Locations));
distanceMetric = distanceMetric.*bJ;
distanceMetric = double(distanceMetric);
distanceMetric(distanceMetric ==0) = Inf;
distanceMetric(1:1+size(distanceMetric,1):end) = 0;
[metric,path] = FloydWarshall(distanceMetric);
metrics(:,:,5) = metric;
paths(:,:,5) = path;

bJ = double(bJ);
bJ(bJ ==0) = Inf;
bJ(1:1+size(bJ,1):end) = 0;
[metric,path] = FloydWarshall(bJ);
metrics(:,:,6) = metric;
paths(:,:,6) = path;

bJ = JTree ~= 0;
distanceMetric = squareform(pdist(Locations));
distanceMetric = distanceMetric.*bJ;
distanceMetric = double(distanceMetric);
distanceMetric(distanceMetric ==0) = Inf;
distanceMetric(1:1+size(distanceMetric,1):end) = 0;
[metric,path] = FloydWarshall(distanceMetric);
metrics(:,:,7) = metric;
paths(:,:,7) = path;

bJ = double(bJ);
bJ(bJ ==0) = Inf;
bJ(1:1+size(bJ,1):end) = 0;
[metric,path] = FloydWarshall(bJ);
metrics(:,:,8) = metric;
paths(:,:,8) = path;

bJ = JGSP ~= 0;
distanceMetric = squareform(pdist(Locations));
distanceMetric = distanceMetric.*bJ;
distanceMetric = double(distanceMetric);
distanceMetric(distanceMetric ==0) = Inf;
distanceMetric(1:1+size(distanceMetric,1):end) = 0;
[metric,path] = FloydWarshall(distanceMetric);
metrics(:,:,9) = metric;
paths(:,:,9) = path;

bJ = double(bJ);
bJ(bJ ==0) = Inf;
bJ(1:1+size(bJ,1):end) = 0;
[metric,path] = FloydWarshall(bJ);
metrics(:,:,10) = metric;
paths(:,:,10) = path;




function [D,P] = FloydWarshall(D)
    %The input weight (or initial distance) matrix must have Inf values where the nodes aren't connected and 0's on the diagonal.
    % Outputs are the shortpaths' distance matrix D, and predecessor's matrix P such that P(i,j) is the node before j on the shortest path from i to j, so if you want to build the paths you have to read P backwards.
    % By Giorgos Dim   
    prevD = D;
    P = zeros(size(D));
    for k = 1:length(D)
        k
        D = min(D,D(:,k) + D(k,:));
        P(D<prevD) = k;
        prevD = D;
    end
end