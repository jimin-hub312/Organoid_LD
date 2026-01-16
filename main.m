%% What this script does (per timepoint):
%   1) Load pre-cropped RI volumes (RI_crop_cell) provided as datasets.
%   2) Generate multi-threshold regional-maxima candidates 
%   3) Filter candidates by 3D morphology + 2D Feret ratio per Z-slice
%   4) Build final LD mask + save mask, MIP visualization, and per-LD stats

% ---------------------- INPUT ----------------------
inMatPath = '/path/to/input';  % MAT containing RI_crop_cell, tpNames, etc.
resultsDir = '/path/to/output/'; % output directory for masks/stats/MIPs
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

% ---------------------- Load MAT ----------------------
% Expected variables in MAT:
%   RI_crop_cell : cell array, each cell is RI_crop (Y x X x Z), double
%   tpNames      : cell array of strings/chars, e.g., {'000000','000001',...}
%   voxelSize    : [vx vy vz] in µm
%   RI_bg        : background RI (e.g., 1.337)
%   alpha        : RI-to-concentration conversion factor
%   tcfPath      : original TCF path (used only to define baseName)

L = load(inMatPath, 'RI_crop_cell', 'tpNames', 'voxelSize', 'RI_bg', 'alpha', 'tcfPath');

RI_crop_cell = L.RI_crop_cell; 
tpNames      = L.tpNames;

voxelSize = L.voxelSize;
RI_bg     = L.RI_bg;
alpha     = L.alpha;
voxelVol  = prod(voxelSize);

% ---------------------- Parameters ----------------------
thresh    = [0.005 0.007 0.01 0.012 0.015 0.017 0.02 0.025 0.03];

% 3D morphology gates
compactThr  = 4.4;   
solidityThr = 0.50; 
aspectThr   = 2.0;   

strength = 0.3;

% base name
[~, baseName, ~] = fileparts(L.tcfPath);

% ---------------------- Loop over timepoints (tp0~tp4) ----------------------
for ti = 1:numel(RI_crop_cell)  
    tpName = tpNames{ti};
    
    fname = sprintf('%s_%s_masks', baseName, tpName);

    maskFile = fullfile(resultsDir, [fname '_LDmask.mat']);
    statFile = fullfile(resultsDir, [fname '_LDstats.mat']);
    mipFile  = fullfile(resultsDir, sprintf('%s_MIP.tiff', fname));

    fprintf('Processing %s (tp=%s)\n', baseName, tpName);

    RI_crop = RI_crop_cell{ti}; 
    if isempty(RI_crop)
        warning('RI_crop is empty at tp=%s. Skipping.', tpName);
        continue;
    end
    if any(isnan(RI_crop(:)))
        warning('RI_crop contains NaN at tp=%s. Skipping.', tpName);
        continue;
    end

    % Candidate detection: Regional maxima for multiple thresholds
    RImax = false([size(RI_crop), numel(thresh)]);
    for i = 1:numel(thresh)
        RImax(:, :, :, i) = imextendedmax(RI_crop, thresh(i));
    end

    % Morphological filtering and LD mask construction 
    LD = zeros(size(RI_crop));
    LD_feret = zeros(size(RI_crop));  

    res = [];

    for idx = length(thresh):-1:1
        cc    = bwconncomp(RImax(:,:,:,idx));
        stats = regionprops3(cc, 'SurfaceArea','Volume','Solidity','BoundingBox','VoxelIdxList','Image');

        for row = 1:size(stats,1)
            % 1) Compactness
            comp = stats.SurfaceArea(row).^(1/2) / stats.Volume(row).^(1/3);
            if comp > compactThr
                continue
            end

            % 2) Solidity
            if stats.Solidity(row) < solidityThr
                continue
            end

            % 3) X–Y Aspect ratio
            box      = stats.BoundingBox(row,:);
            box_size = box(4:6);
            arXY = max(box_size(1),box_size(2)) / min(box_size(1),box_size(2));
            if arXY > aspectThr
                continue
            end

            LDimg = stats.Image{row};

            isferetsuccessful = true;
            for z = 1:size(LDimg,3)
                LDimg2D = LDimg(:,:,z);
                stats2d = regionprops(LDimg2D, 'MinFeretProperties','MaxFeretProperties');

                maxD = [stats2d.MaxFeretDiameter];
                minD = [stats2d.MinFeretDiameter];

                feretRatio = maxD(1) / (minD(1)+2);

                if feretRatio > 1.5
                    isferetsuccessful = false;
                    continue
                end
            end

            if ~isferetsuccessful
                continue
            end
           
            res = [res;thresh(idx) stats.Volume(row) comp stats.Solidity(row)]; %#ok<AGROW>

            LD(stats.VoxelIdxList{row}) = max(thresh(idx), LD(stats.VoxelIdxList{row}));
        end
    end

    LD = imfill(LD, 'holes');
    LDmask = LD > 0;

    % Post-filter: Filter small vs large objects, produce finalMask
    cc  = bwconncomp(LDmask);
    finalMask = false(size(LDmask));
    for j = 1:cc.NumObjects
        vox = cc.PixelIdxList{j};
        tmp = false(size(LDmask)); tmp(vox) = true;
        mip2D = max(tmp, [], 3);
        rp = regionprops(mip2D, 'MajorAxisLength');
        if isempty(rp), continue; end
        if rp.MajorAxisLength <= 5
            keep = vox(LD(vox) >= 0.01);
        else
            keep = vox;
        end
        finalMask(keep) = true;
    end
    LDmask_filtered = finalMask;

    % Save mask
    save(maskFile, 'LDmask_filtered');

    % Visualization MIP 
    highlight = LDmask_filtered;

    volNorm = mat2gray(RI_crop);
    RGB     = repmat(volNorm, [1 1 1 3]);
    RGB(:,:,:,2) = max(RGB(:,:,:,2) - strength*highlight, 0);
    RGB(:,:,:,3) = max(RGB(:,:,:,3) - strength*highlight, 0);

    RGB_high = squeeze(max(RGB, [], 3));

    RGB_high = im2uint8(RGB_high);
    mipGray  = im2uint8(repmat(mat2gray(squeeze(max(RI_crop,[],3))), [1 1 3]));

    merged   = [mipGray, RGB_high];

    outName = string(fname) + "_MIP.tiff";
    outPath = fullfile(resultsDir, outName);
    imwrite(merged, outPath);

    % Compute LD statistics
    [lbl, N] = bwlabeln(LDmask_filtered);
    stats = regionprops3(lbl, RI_crop, 'VoxelIdxList', 'Volume', 'SurfaceArea', ...
        'MeanIntensity', 'MinIntensity', 'MaxIntensity');

    LD_Index   = (1:N)';
    Conc       = nan(N,1);
    Drymass    = nan(N,1);
    ESD        = nan(N,1);
    Sphericity = nan(N,1);

    SurfaceArea_um2 = stats.SurfaceArea * prod(voxelSize(1:2));
    Volume_um3      = stats.Volume     * voxelVol;

    RI_max  = stats.MaxIntensity;
    RI_mean = stats.MeanIntensity;
    RI_min  = stats.MinIntensity;

    RI_std = nan(N,1);
    for j = 1:N
        vals = RI_crop(stats.VoxelIdxList{j}); RI_std(j) = std(vals);
        if RI_mean(j) <= RI_bg, continue; end
        V = Volume_um3(j);
        c = (RI_mean(j)-RI_bg)/alpha;
        if c <= 0, continue; end
        Conc(j)       = c;
        Drymass(j)    = c * V;
        ESD(j)        = (6*V/pi)^(1/3);
        Sphericity(j) = (pi^(1/3)) * (6*V)^(2/3) / SurfaceArea_um2(j);
    end

    varNames = {'LD_Index','Concentration','Drymass','ESD','RI_max','RI_mean','RI_min','RI_std', ...
        'Sphericity','SurfaceArea_um2','Volume_um3'};
    Results = table(LD_Index, Conc, Drymass, ESD, RI_max, RI_mean, RI_min, RI_std, ...
        Sphericity, SurfaceArea_um2, Volume_um3, 'VariableNames', varNames);
    Results(isnan(Results.Concentration), :) = [];

    save(statFile, 'Results');

    fprintf('Completed tp=%s\n', tpName);
end

disp('Done.');
toc;

