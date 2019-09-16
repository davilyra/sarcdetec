clear all
close all
clc

% Initializing general elements
hist_centers=-90:3:90;

% Select the first file of the sequence
disp('Select the first file from the data sequence')
[file,path]=uigetfile({'*.mat';'*.MAT'},'Select MAT File...');
filename = [path file];
disp(' ')

% Defining the number of input datasets
n = input('What is the number of datasets that you wish to combine?\n');

% Preallocating memory (speed up step)
images = cell(n);
all_images = [];
fname = filename(1:end-29);

ii2 = 0;
for ii = str2double(filename(end-28)): str2double(filename(end-28)) + n - 1
    ii2 = ii2 + 1;
    load([fname,num2str(ii),'.tif.sarcDetect.Settings.mat'])
    images{ii2} = rad2deg(nonzero_orientation);

    subplot(ceil(n/3),3,ii2)
    image_hist = hist(images{ii2},hist_centers);
    bar(hist_centers,100*image_hist/sum(image_hist))
    xlabel('Sarcomere Orientation')
    ylabel('Normalized Occurence')
    title(['Image ',num2str(ii)])
    axis([-90 90 0 5])
    
    all_images = [all_images; images{ii2}];
end

all_images_hist = hist(all_images,hist_centers);
all_images_hist = 100*all_images_hist/sum(all_images_hist);
figure
bar(hist_centers, all_images_hist, 1, 'k')
xlabel('Sarcomere Orientation')
ylabel('Normalized Occurence')
title('All Images - Un-Shifted to Mode')
axis([-90 90 0 5])


[~,index] = max(all_images_hist)
shift = hist_centers(index)
all_images_shifted = all_images - shift;
for i = 1:size(all_images_shifted)
    if all_images_shifted(i) < -90
        all_images_shifted(i) = all_images_shifted(i) + 180;
    end
    if all_images_shifted(i) > 90
        all_images_shifted(i) = all_images_shifted(i) - 180;
    end
end

figure
all_images_shifted_hist = hist(all_images_shifted,hist_centers);
all_images_shifted_hist = 100*all_images_shifted_hist/sum(all_images_shifted_hist);
bar(hist_centers,all_images_shifted_hist,1,'k')
xlabel('Sarcomere Orientation')
ylabel('Normalized Occurence')
title('All Images - Shifted to Mode')
axis([-90 90 0 5])

% Saving Data
save all_images_sarc_alignment_unshifted 'all_images' 'hist_centers'
clear all_images
all_images = all_images_shifted;
% save all_images_sarc_alignment_shifted_to_mode 'all_images_shifted'
save all_images_sarc_alignment_shifted_to_mode 'all_images' 'hist_centers'
