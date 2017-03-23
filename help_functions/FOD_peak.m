%% function to find peaks of a spherical function

function [peak_value, peak_idx, peak_theta, peak_phi, peak_pos, peak_pos_final] = FOD_peak(fod, dis, nbhd, thresh, pos, theta, phi)

%%%%% Output
%% peak_value: function value of the peaks
%% peak_idx: index of the peaks in the vector of function values
%% peak_theta,peak_phi: peak locations in angles
%% peak_pos: peak location in cartesian coordinates
%% peak_pos_final: final positions of peaks after grouping close peaks and averaging locations

%%%%% Input
%% fod: spherical function
%% dis: pairwise shperical distance of grid points used in finding nbhd
%% nbhd: number of nbhd grid points of each grid point (including itself: 40 for 2562 equal angle grid points)
%% general rull is to make each neighborhood span around 30 degree
%% thresh: to define cutoff threshhold, thresh*max(fod) is the cutoff value (less than it=disgard)
%% pos, (theta, phi): cartesian (spherical) coordinates of the grid points
    peak_idx = [];
    
    n = size(pos, 2);
    
    % find local peaks in each neighbor hood 
    for i=1:n
        [~, AIdx] = sort(dis(i,:));
%         smallestNElements = ASorted(2:kmin);
        smallestNIdx = AIdx(2:nbhd);
        if(fod(i)>= quantile(fod(smallestNIdx),1))
            peak_idx = [peak_idx, i];
        end     
    end
    
  
    cutoff = thresh*max(fod);
    peak_idx = peak_idx(fod(peak_idx)>cutoff);
    
    peak_pos = pos(:,peak_idx);
    
    idx_dele = [];
    for i=1:size(peak_idx,2)
        for j=(i+1):size(peak_idx,2)
            if (-1e-2)<(peak_pos(:,i)'*peak_pos(:,j)+1) && (peak_pos(:,i)'*peak_pos(:,j)+1)<(1e-2)
                idx_dele = [idx_dele, j];
            end
        end
    end
    
    peak_idx(idx_dele) = [];
    
    peak_pos = pos(:,peak_idx);
    peak_value = fod(peak_idx);
    peak_theta = theta(peak_idx);
    peak_phi = phi(peak_idx);
    
    k = numel(peak_idx);
    pos_cross_norm = zeros(k,k);
    for i=1:k
        for j =1:k
            pos_cross_norm(i,j) = norm(cross(peak_pos(:,i),peak_pos(:,j)));
        end
    end
    aaa = atan2(pos_cross_norm,peak_pos'*peak_pos);
    aaa_dup = aaa*180/pi<5 & aaa*180/pi>-5;
    
    peak_pos_final = [];
    idx_used = [];
    
    for i=1:size(peak_idx,2)
        group = [];
        for j=i:size(peak_idx,2)
            if (~(any(j==idx_used))) && (aaa_dup(i,j)==1)
                group = [group,j];
                idx_used = [idx_used, j];
            end
        end 
        
        if isempty(group) ~= 1
            deno = sum(peak_value(group));

            tmat = repmat(peak_value(group), 3, 1 )./deno;

            peak_pos_temp = sum(peak_pos(:,group).*tmat,2);
            peak_pos_temp = peak_pos_temp./(sum(peak_pos_temp.^2));
            peak_pos_final = [peak_pos_final, peak_pos_temp];
        end
    end
    
end