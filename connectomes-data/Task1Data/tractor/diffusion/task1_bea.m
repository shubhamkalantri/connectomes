% info = load_nii('T1w_vol1.nii');
% parcellation = load_nii('parcellation.nii');
% view_nii(info)
% cd ../../connectomes/

dti = load_nii('dti_FA.nii');
threshold = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8];
% Initialize a cell array to store thresholded images
FA_list = cell(1, 8);
for i = 1:8
    thres_val = threshold(i);  % this is now a numeric value
    FA = dti.img;              % get the FA image data
    % Create a thresholded image: set voxels below the threshold to zero
    mask = FA;
    mask(mask < thres_val) = 0;
    FA_list{i} = mask;
end

% Visualize the middle slice of each thresholded image in subplots
figure;
numThresholds = length(threshold);
for i = 1:numThresholds
    subplot(2,4,i);
    % For a 3D volume, take the middle slice along the third dimension
    middleSlice = round(size(FA_list{i}, 3) / 2);
    imshow(FA_list{i}(:,:,middleSlice), []);
    title(sprintf('Threshold = %.1f', threshold(i)));
end
%% 1.3 bct to calculate edges
%  Load the structural connectomes in the connectomes directory and use bct
% edge density, mean shortest path,efficiency and mean clustering coefficient,
% cd ../../connectomes/
fa_tables = cell(1,8);
edge_densities = zeros(1, 8);
mspath = zeros(1,8); 
mefficiency = zeros(1,8);
cluster_coef = zeros(1,8);

for i=1:8
    filename = sprintf('FA.%d_graph.csv', i);
    fa = csvread(filename,3,0);
    fa_tables{i} = fa;
    % maybe use clique to find biggest componen
    G = graph(fa); % Convert to graph object
    comp = conncomp(G);
    largestComp = mode(comp);  
    nodesToKeep = find(comp == largestComp); % Extract nodes in the largest component
    maxMatrix = fa(nodesToKeep, nodesToKeep); 
    edge_densities(i) = density_und(fa);
    % get distance matrix
    D = distance_bin(fa);
    [mspath(i),mefficiency(i)] =charpath(D);
    cluster_coef(i) = mean(clustering_coef_bu(fa));
end
% not finding the max clique gives this 
%  ed 0.1111    0.1093    0.0957    0.0904    0.0812    0.0663    0.0443    0.0281
%  mspath  2.6883    2.7884       Inf       Inf       Inf       Inf       Inf       Inf
%  meffic  0.4413    0.4322    0.3894    0.3845    0.3430    0.2904    0.1825    0.0756
%  cluster coef all zeros 00.4880    0.5143    0.4517    0.4993    0.3980    0.3282    0.2549    0.1696

% after using conncomp and get max component 
% ed           0.1111    0.1093    0.1012    0.0932    0.0858    0.0825    0.0833    0.1225
% mspath       2.6883    2.7884    2.9478    3.1465    3.5823    3.5885    3.7491    3.1897
% meffic       0.4413      0.4322    0.4131    0.3962    0.3638    0.3615    0.3509    0.4041
% cluster coef 0.4880    0.5143    0.4654    0.5067    0.4101    0.3658    0.3538    0.2710










