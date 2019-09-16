clear all
close all

load aact1.tif.sarcDetect.Settings.mat
image1=rad2deg(nonzero_orientation);
load aact2.tif.sarcDetect.Settings.mat
image2=rad2deg(nonzero_orientation);
load aact3.tif.sarcDetect.Settings.mat
image3=rad2deg(nonzero_orientation);
load aact4.tif.sarcDetect.Settings.mat
image4=rad2deg(nonzero_orientation);
load aact5.tif.sarcDetect.Settings.mat
image5=rad2deg(nonzero_orientation);

hist_centers=-90:3:90;

subplot(2,3,1)
image1_hist=hist(image1,hist_centers);
bar(hist_centers,100*image1_hist/sum(image1_hist))
xlabel('Sarcomere Orientation')
ylabel('Normalized Occurence')
title('Image1')
axis([-90 90 0 5])
subplot(2,3,2)
image2_hist=hist(image2,hist_centers);
bar(hist_centers,100*image2_hist/sum(image2_hist))
xlabel('Sarcomere Orientation')
ylabel('Normalized Occurence')
title('Image2')
axis([-90 90 0 5])
subplot(2,3,3)
image3_hist=hist(image3,hist_centers);
bar(hist_centers,100*image3_hist/sum(image3_hist))
xlabel('Sarcomere Orientation')
ylabel('Normalized Occurence')
title('Image3')
axis([-90 90 0 5])
subplot(2,3,4)
image4_hist=hist(image4,hist_centers);
bar(hist_centers,100*image4_hist/sum(image4_hist))
xlabel('Sarcomere Orientation')
ylabel('Normalized Occurence')
title('Image4')
axis([-90 90 0 5])
subplot(2,3,5)
image5_hist=hist(image5,hist_centers);
bar(hist_centers,100*image5_hist/sum(image5_hist))
xlabel('Sarcomere Orientation')
ylabel('Normalized Occurence')
title('Image5')
axis([-90 90 0 5])

all_images=[image1;image2;image3;image4;image5];
all_images_hist=hist(all_images,hist_centers);
all_images_hist=all_images_hist/sum(all_images_hist);
bar(hist_centers,all_images_hist)

[mode,index]=max(all_images_hist)
shift=hist_centers(index)
all_images_shifted=all_images-shift;
for i=1:size(all_images_shifted)
    if all_images_shifted(i)<-90
        all_images_shifted(i)=all_images_shifted(i)+180;
    end
    if all_images_shifted(i)>90
        all_images_shifted(i)=all_images_shifted(i)-180;
    end
end

all_images_hist=hist(all_images,hist_centers);
all_images_hist=100*all_images_hist/sum(all_images_hist);
figure
bar(hist_centers,all_images_hist,1,'k')
xlabel('Sarcomere Orientation')
ylabel('Normalized Occurence')
title('All Images - Un-Shifted to Mode')
axis([-90 90 0 5])
figure
all_images_shifted_hist=hist(all_images_shifted,hist_centers);
all_images_shifted_hist=100*all_images_shifted_hist/sum(all_images_shifted_hist);
bar(hist_centers,all_images_shifted_hist,1,'k')
xlabel('Sarcomere Orientation')
ylabel('Normalized Occurence')
title('All Images - Shifted to Mode')
axis([-90 90 0 5])
save all_images_sarc_alignment_unshifted 'all_images'
save all_images_sarc_alignment_shifted_to_mode 'all_images_shifted'