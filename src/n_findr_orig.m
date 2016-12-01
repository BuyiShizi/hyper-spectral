%% N-FINDR function for endmember extraction
%-----------------------------------------------------------------------------------
    % this function is a realization for a new endmember extraction method, which is
    % similiar to the original N-FINDR method.
%-----------------------------------------------------------------------------------
function [endmember_index, V, M] = n_findr_orig (hyspectral_data, endmember_number)
    % hyperspectral_data:   input hyperspectral_data
    % endmember_number:     endmember endmember_number
    % endmember_index:      output index of endmember in original hyperspectral_data
    % V:                    output volumn conposed by endmember
    % M:                    media data
    [row_hyper, col_hyper] = size (hyspectral_data);
    A = zeros (col_hyper, endmember_number); % store diff vector
    V_max = 0; 
    flag_skip = 0;
    endmember_idx = zeros(1, endmember_number);
    endmember_vector = zeros (col_hyper, endmember_number); % store endmember_vector in every iteration
    iter_1 = 1;
    iter_2 = 1;
    iter_3 = 1;
    M = -1 * ones (endmember_number, row_hyper);
    m_i = 0;
    [man_val, min_indx]  = min( hyspectral_data(:,1));
    endmember_idx(1) = min_indx;
    endmember_vector(:,1) = hyspectral_data(min_indx, :)';
    for i = 2:endmember_number % iteration for endmember
        A = zeros (col_hyper,i);
        if (2==i) % the first endmember extraction
            m_i = 0;
            for j = 1:row_hyper % traverse all hyperspectral data
                A = hyspectral_data(j,1:i-1)';
                AAT = det (A' * A);
                AAT = sqrt (AAT); % volumn of endmember
                if AAT > V_max % update endmember
                    V_max = AAT;
                    endmember_vector(:,i) = hyspectral_data(j, :)';
                    endmember_idx(i) = j;
                    endmember_1_change(iter_1) = j;
                    iter_1 = iter_1 + 1;
                    
                    m_i = m_i + 1;
                    M (i,m_i) = j;
                end
            end
        else
            m_i = 0;
            for ii = 1:i-1
                    A(:,ii) = endmember_vector(:,ii);
%                     A(:,ii) = hyspectral_data(j,:)' - endmember_vector(:,ii);
            end
            for j = 1:row_hyper
                for jj = 1:i 
                    if (j==endmember_idx(jj)) || sum(hyspectral_data(j,:)'-endmember_vector(:,jj))==0
                        flag_skip = 1;
                        break;
                    end
                end
                if flag_skip==1
                    flag_skip = 0;
                    continue;
                end
                A(:,i) = hyspectral_data(j,:)';
                AAT = det (A' * A);
                AAT = sqrt (AAT) / factorial(i-2);
                if AAT > V_max
                    V_max = AAT;
                    endmember_vector(:,i) = hyspectral_data(j,:)';
                    endmember_idx(i) = j;
                    if i==2
                        endmember_3_change(iter_3) = j;
                        iter_3 = iter_3 + 1;
                        
                        m_i = m_i + 1;
                        M (i,m_i) = j;
                    else
                        endmember_3_change(iter_3) = j;
                        iter_3 = iter_3 + 1;
                        
                        m_i = m_i + 1;
                        M (i,m_i) = j;
                    end
                end
            end
        end
        V_max = 0;
    end
    endmember_index = endmember_idx;
    V = V_max;
    
end
