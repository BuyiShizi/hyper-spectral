function [loc, V_max] = finder(pca_img, nlength, Nend)
%% N_FINDR algorigthm

% randomly select pixels as initial endmembers.
ind = unidrnd(nlength, 1, Nend);%the process begin with a random set of pixel as endmembers 1 * Nend
E   = ones(Nend);% Nend * Nend
E(2:end,:) = pca_img(ind,:)';% E stores iteration endmembers
dentor = factorial(Nend-1);%
V_max  = abs(det(E))/dentor;% % calculate n-dementional volumn

% find the largest volume, and set the corresponding pixels as
% endmembers.
for i = 1:Nend  % loop over each endmember.
    i_max = ind(i);
    for j = 1:nlength  % loop over each pixel.
        E(2:end,i) = pca_img(j,:)';
        V = abs(det(E))/dentor;
        k = V>V_max;
        i_max(k) = j;
        V_max(k) = V;
    end
    E(2:end,i) = pca_img(i_max,:)'; % see <N_FINDR> Page 3
    ind(i) = i_max;
end
loc = ind;

% end function finder().
end