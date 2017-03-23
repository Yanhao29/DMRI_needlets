%% function that changes index of real data to fit matrix tradition
%% in real dti data, dim1 is from left to right (column-wise), dim2 is from bot to top (row-wise)
function dti_data_reindex = dti_data_reindex(data)
    
    data = squeeze(data);
    n1 = size(data,1);
    n2 = size(data,2);
    if(length(size(data))==3)
        n3 = size(data,3); 
        data_temp = zeros(n2,n1,n3);
        for k1 = 1:n2
            for k2 = 1:n1
                data_temp(k1,k2,:) = data(k2,n2+1-k1,:);
            end
        end
    else
        data_temp = zeros(n2,n1);
        for k1 = 1:n2
            for k2 = 1:n1
                data_temp(k1,k2) = data(k2,n2+1-k1);
            end
        end
    end
    
    for k1 = 1:n2
        for k2 = 1:n1
            data_temp(k1,k2,:) = data(k2,n2+1-k1,:);
        end
    end
    dti_data_reindex = data_temp;  

end