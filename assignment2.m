%cmpt412 Assignment2
%author: Janet Xuan
%
% image=imread('Cig01.JPG'); % 
% image=imread('Cig03.JPG'); %area too small
% image=imread('Cig05.JPG'); % diff cig color, too dark edge not seen
% image=imread('Cig07.JPG'); %noise
% image=imread('Cig08.JPG'); %too dark, edge not seen
% image=imread('Cig09.JPG'); 
% image=imread('Cig13.JPG'); %noise
image=imread('Cig_on_Orange1.JPG'); %


bw=(double(image(:,:,1))+double(image(:,:,2))+double(image(:,:,3)))/(3*255);
%%normalized to 0-1 range


% % https://www.mathworks.com/help/images/examples/detecting-a-cell-using-image-segmentation.html
figure, imshow(bw)
[~, threshold] = edge(bw, 'Prewitt');
fudgeFactor =1;
BWs = edge(bw,'Prewitt', threshold * fudgeFactor);
figure, imshow(BWs), title('binary gradient mask');

se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);

BWsdil = imdilate(BWs, [se90 se0]);
figure, imshow(BWsdil), title('dilated gradient mask');

BWdfill = imfill(BWsdil, 'holes');
figure, imshow(BWdfill);
title('binary image with filled holes');

BWnobord = imclearborder(BWdfill, 4);
figure, imshow(BWnobord), title('cleared border image');

seD = strel('diamond',1);
BWfinal = imerode(BWnobord,seD);
BWfinal = imerode(BWfinal,seD);
figure, imshow(BWfinal), title('segmented image');


%%%%%%%%%       locate    %%%%%%%

BW3 = bwareaopen(BWfinal, 50);

% 
BWlabel = bwlabel(BW3); %label connected objects
STATS = regionprops(BWlabel, 'Area' );

cig_areas = [STATS.Area];

BW4 = bwareaopen(BW3, round(mean(cig_areas)/2));


% cc = bwconncomp(BWlabel, 4);

% cig = false(size(BWlabel));


% for k = 1 : length(STATS)
%     idx = STATS(k).PixelIdxList;
%     if STATS(k).Area > mean(cig_areas)
%         cig(STATS.PixelIdxList{idx}) = true;
%     end
%     
%         
%     % Now do something with them.....
% end


cigs = regionprops(BW4,'centroid'); %find center of each cigs
centroids = cat(1, cigs.Centroid); %concatenate the center coordinates
figure, imshow(BW4); hold on;

plot(centroids(:,1),centroids(:,2),'go', 'MarkerFaceColor', 'g','MarkerSize', 10)
disp('Cigs Center Coordinates: ')
disp(centroids)









