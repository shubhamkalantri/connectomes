function [structural,functional] = LoadConnectomes(path)
%LOADCONNECTOMES Read in structural and functional connectomes.
    data = zeros(2, 19, 68, 68);
    groupnames = ["WFA", "rsfMRI"];
    for connectomeType = 1:numel(groupnames)
        group = groupnames(connectomeType);
        files = dir(sprintf("%s/*_%s_68.csv", path, group));
        assert(numel(files) == 19, "b ruh");
        [~,ind] = sort({files.name});
        files = files(ind);
        % NOTE: subjects are numbered 32 to 50, who's alphabetical ordering is
        %       the same as the numerical ordering, so subjects are indeed
        %       ordered by their ID number.
        for subject = 1:19
            filename = fullfile(path, files(subject).name);
            data(connectomeType, subject, :, :) = readmatrix(filename);
        end
    end
    structural = squeeze(data(1, :, :, :));
    functional = squeeze(data(2, :, :, :));
end