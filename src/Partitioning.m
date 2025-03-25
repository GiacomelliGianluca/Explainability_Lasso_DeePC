function [label, Centr] = Partitioning(Data, thr)
% Partitioning of Data based on the threshold thr

% label: label of the data (2: above threshold, 1: under threshold)
% Centr: Centroids of the data

% Allocating label vector
label = NaN*ones(length(Data),1);
% Allocating Centroid vector 
Centr = NaN*ones(2,1);

for i = 1 : length(Data)
    if Data(i) >= thr
        label(i) = 2;
    else
        label(i) = 1;
    end
end

for i = 1 : 2
    Centr(i) = mean(Data(label==i));
end



end