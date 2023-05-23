clc;clear;close all

%set the load and save path
path  = 'D:\Individualinform\YiCheng\FeedbackforWH\image2';
dirname = "CPTAC-GBM";
raw_loc= fullfile(path, "clear_data_lower", dirname);
seg_loc= fullfile(path, "out", dirname);
save_loc = 'video';

FL = dir(fullfile(seg_loc, '*.nii.gz'));

% Set color and the transparent
colorMapWT = [0.4 0.2 0.8];
colorMapTC = [0.2, 0.4, 0.8];
colorMapET = [0.8, 0.4, 0.2];
colorMaps = [colorMapWT;
             colorMapTC;
             colorMapET];
alpha = 0.3;

% make directory
if ~exist('video', 'dir')
    mkdir(save_loc)
end

for k = 1
    filename = FL(k).name;
    prefix = replace(filename, '.nii.gz', '');   
    img_path = fullfile(raw_loc, prefix, [prefix '_t1ce.nii.gz']);
    seg_path = fullfile(seg_loc, filename);
    save_path = fullfile(save_loc, prefix);
    Out = Write(img_path, seg_path, save_path);
end

%% Function 
function Out = Write(raw_path, seg_path, save_path, colorMaps, alpha)

if ~exist('alpha', 'var')
    alpha = 0.3;
end

if ~exist('colorMaps', 'var')
    colorMapWT = [0.4 0.2 0.8];
    colorMapTC = [0.2, 0.4, 0.8];
    colorMapET = [0.8, 0.4, 0.2];
    colorMaps = [colorMapWT;
                 colorMapTC;
                 colorMapET];
end

img = niftiread(raw_path);
seg = double(niftiread(seg_path));
out = MakeFrame(img, seg, colorMaps, alpha);
Out = [out.newImg, out.WTImg, out.TCImg, out.ETImg;];
WriteVideo(Out, [save_path '.avi']);
end


function [Out] = MakeFrame(img, pred, colorMaps, alpha)
    colorMapWT = colorMaps(1,:);
    colorMapTC = colorMaps(2,:);
    colorMapET = colorMaps(3,:);
    [newImg] = MinMaxNormal(img);
    [Out.newImg] = AddRGBChannel(newImg);
    [Out.WTImg] = Overlay(Out.newImg, ismember(pred, [1,2,4]), colorMapWT, alpha);
    [Out.TCImg] = Overlay(Out.newImg, ismember(pred, [1,4]), colorMapTC, alpha);
    [Out.ETImg] = Overlay(Out.newImg, ismember(pred, [4]), colorMapET, alpha);
end

function [newImg] = MinMaxNormal(image)
m = min(image, [], 'all');
M = max(image, [], 'all');
newImg = (image - m) / (M-m);

end

function [newImg] = AddRGBChannel(image)
newImg = repmat(image, 1, 1, 1, 3);

end

function [outImg] = Overlay(downImg, upImg, colorMap, alpha)
mask = downImg;
nonZero = AddRGBChannel(upImg ~= 0);
imgSize = size(downImg);
colorMask = zeros(imgSize);
colorMask(:, :, :, 1) = colorMap(1);
colorMask(:, :, :, 2) = colorMap(2);
colorMask(:, :, :, 3) = colorMap(3);
colorMask =  colorMask .* nonZero;
mask(nonZero) = colorMask(nonZero);
outImg = alpha*mask + (1-alpha)*downImg;

end

function WriteVideo(image, videoName)

v = VideoWriter([videoName]);
[~, ~, nz, ~] = size(image);

open(v);

for i= 1:nz
    slice = squeeze(image(:, :, i, :));
    writeVideo(v, slice);
end

close(v);

end