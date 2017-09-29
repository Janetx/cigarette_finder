%cmpt412 Assignment2
%author: Janet Xuan
%
% image=imread('Cig01.JPG');
% image=imread('Cig03.JPG');
% image=imread('Cig05.JPG');
% image=imread('Cig07.JPG');
% image=imread('Cig08.JPG');
% image=imread('Cig09.JPG');
image=imread('Cig13.JPG');
% image=imread('Cig_on_Orange1.JPG');

bw=(double(image(:,:,1))+double(image(:,:,2))+double(image(:,:,3)))/(3*255);
%%normalized to 0-1 range

%%%%*********   Simple gradient-based edge detection   *****************
Threshold = 0.1;
cx=conv2(bw,[-1 1],'same');  %Horizontal derivative
cy=conv2(bw,[-1 1]','same'); %Vertical derivative
magsq=cx.^2+cy.^2;   %Magnitude of the gradient
edgesdetected=zeros(size(magsq));
edgesdetected(magsq/max(magsq(:))>Threshold)=1;  %Threshold on gradient
imshow(edgesdetected) %%% good picture

% **************    locate     *****************
 BW = bwareaopen(edgesdetected, 50);
% 
BWlabel = bwlabel(edgesdetected); %label connected objects
cigs = regionprops(BWlabel,'centroid'); %find center of each pencil
centroids = cat(1, cigs.Centroid); %concatenate the center coordinates
figure, imshow(BW); hold on;

plot(centroids(:,1),centroids(:,2),'go', 'MarkerFaceColor', 'g','MarkerSize', 10)
disp('Cigs Center Coordinates: ')
disp(centroids)

% input('Simple gradient edge detection result. Hit any key to continue')

%Run examples of Canny edge detection
% ec = edge(bw,'canny');
% figure; imshow(ec);
% ec = edge(bw,'canny', .1, 1);
% figure; imshow(ec);
% ec = edge(bw,'canny', .1, 10);
% figure; imshow(ec);

%Examples of zero-crossing of Laplacian edge detection

% ez = edge(bw,'zerocross');
% 
%     eg = edge(bw,'log'); 
%     figure; imshow(eg);
    
%     input('Canny edge result. Hit any key to continue')
    
[eg thresh]= edge(bw,'log');


[eg thresh]= edge(bw,'log', thresh, 3);
figure; imshow(eg);
% [eg thresh]= edge(bw,'log', thresh, 4);
% figure; imshow(eg);
% [eg thresh]= edge(bw,'log', thresh, 5);
% figure; imshow(eg);
disp('LoG results for various thresholds')


% % % **********************    filtering      *****************************

figure; imshow(bw);
%Uniform averaging mask
% h=fspecial('average',15); %h = fspecial('average', hsize) returns
                          % an averaging filter h of size hsize.
% bg=filter2(h, bw); %Like conv2 but correlation (i.e., filter not reflected)

% figure; imshow(bg);
% title('Filtered by Uniform Mask')
% figure; surf(h);
% title('Uniform mask')
% input('Hit a key to continue');

% h = fspecial('gaussian',13,3);
% bg=filter2(h, bw);
% figure; imshow(bg);
% title('Filtered by Gaussian Mask with Sigma = 3')
% figure; surf(h);
% title('Gaussian Mask Sigma = 3')
% input('Hit a key to continue');

% h = fspecial('gaussian',31,7);
% bg=filter2(h, bw);
% figure; imshow(bg);
% title('Filtered by Gaussian Mask with Sigma = 7')

%Display a plot of the convolution mask h
% figure; surf(h);
% input('Hit a key to continue');
%Examples of the laplacian of Gaussian
%Sigma of 3. Masksize of 15
h = fspecial('log',15,3);
bg=filter2(h, bw);
% figure; surf(h);  %Display the mask
% title('Laplacian of Gaussian Mask (Sigma = 15)')
% input('Hit a key to continue');
%Can't display the result of the convolution directly because it has
%negative values.
%Add an offset and then scale to 0-1 range for display purposes
bgs=bg+.5;
bgs=bgs-min(bgs(:));
bgs=bgs/max(bgs(:));
figure; imshow(bgs);
title('Filtered by Laplacian of Gaussian Mask (plus offset) Sigma = 15')
% input('Hit a key to continue');

% h = fspecial('log',31,5);  %Sigma 5, mask size of 31
% bg=filter2(h, bw);
% figure; surf(h);
% title('Laplacian of Gaussian Mask (Sigma = 5)')
% input('Hit a key to continue');
% bgs=bg+.5;
% bgs=bgs-min(bgs(:));
% bgs=bgs/max(bgs(:));
% figure; imshow(bgs);
% title('Filtered by Laplacian of Gaussian Mask (plus offset) Sigma = 5')


%
% % *************   edge detection after filtering   *****************

% Threshold = 0.1;
% cx=conv2(bgs,[-1 1],'same');  %Horizontal derivative
% cy=conv2(bgs,[-1 1]','same'); %Vertical derivative
% magsq=cx.^2+cy.^2;   %Magnitude of the gradient
% edgesdetected=zeros(size(magsq));
% edgesdetected(magsq/max(magsq(:))>Threshold)=1;  %Threshold on gradient
% figure; imshow(edgesdetected) %%% good picture


% % % %  *****************************   grayscale method %%%%%%%%%%%
%
% 
% figure;
% imshow(image);
% 
% bw=(double(image(:,:,1))+double(image(:,:,2))+double(image(:,:,3)))/(3*255);
% figure;  %open new figure window
% imshow(bw);
% 
% % https://www.mathworks.com/help/images/structuring-elements.html
% background = imopen(bw,strel('diamond',10));
% greyscale2 = bw - background; %remove the background from greyscale
% figure;  %open new figure window
% imshow(greyscale2);
% 
% greyscale3 = imadjust(greyscale2, [0.3 0.4]);
% figure, imshow(greyscale3), title('Increased Contrast');
% 
% level = graythresh(greyscale3);
% BW = im2bw(greyscale3,level);
% figure, imshow(BW);
% 
% BW = bwareaopen(BW, 50);
% 
% BWlabel = bwlabel(BW); %label connected objects
% pencils = regionprops(BWlabel,'centroid'); %find center of each pencil
% centroids = cat(1, pencils.Centroid); %concatenate the center coordinates
% figure, imshow(BW); hold on;
% 
% plot(centroids(:,1),centroids(:,2),'go', 'MarkerFaceColor', 'g','MarkerSize', 10)
% disp('Pencil Center Coordinates: ')
% disp(centroids)








% dt = 110/255;  %Set threshold for being dark. Chosen experimentally
% dark=find(bw<dt); 
% 
% tw=-ones(size(bw));   %Make up a new array full of -1's that's the same
%                       %dimensions as input image bw.
% tw(dark)=1;      %Set all pixels of dark (pennies we hope) indices to +1
% 
% figure; imshow(tw); %Two or more commands can go on one line






