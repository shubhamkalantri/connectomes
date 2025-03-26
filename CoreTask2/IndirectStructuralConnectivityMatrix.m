function T = IndirectStructuralConnectivityMatrix(structural)
%INDIRECTSTRUCTURALCONNECTIVITYMATRIX Calculates the indirect structural
% connectivity matrix given a structural connectome.
    arguments
        structural (:, :, :) double
    end
    T = zeros(size(structural));
    num_subjects = size(structural, 1);
    num_vertices = size(structural, 2);
    for subject = 1:num_subjects
        for i = 1:num_vertices
            for j = 1:num_vertices
                maxVal = -Inf;
                for k = 1:num_vertices
                    if k == i || k == j
                        % skipping self edges
                        continue
                    end
                    s_ik = structural(subject, i, k);
                    s_jk = structural(subject, j, k);
                    % s_ik, s_jk != 0
                    if s_ik == 0 || s_jk == 0
                        continue
                    end
                    % "greatest minimum weight"
                    minimum = min(s_ik, s_jk);
                    maxVal = max(minimum, maxVal);
                end
                if maxVal == -Inf
                    % there are no two-step chains
                    continue
                end
                T(subject, i, j) = maxVal;
            end
        end
    end
end

