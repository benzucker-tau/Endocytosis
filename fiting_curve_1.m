% fiting_curve_1.m

%% using new image analysis tricks

folder      = 'D:\Ben\EVOLVER\chromaffin endocytosis\video' ;
filename    = 'media3.mp4';
%filename    = '03avi .mov'; 

v = VideoReader(fullfile(folder,filename));
numofframes = 8;
totalframes = floor(v.Duration*v.FrameRate);

counter = 0;
counter2 = 0;
% snap = [];
video = readFrame(v);  %% get rids out of the first frame ( should consider running this on a duplicate)
snap = zeros(numofframes,size(video,1),size(video,2),size(video,3)); % ,size(video,3));
while hasFrame(v)
    counter=counter+1 ;
    video = readFrame(v);
    if rem(counter,floor((totalframes)/(numofframes-1))) == 1
        counter2 = counter2+1;
     %   onecolor = video(:,:,2);
      %  snap(counter2,:,:) = onecolor;
        snap(counter2,:,:,:) = video ; %rgb2gray(video) ;
        a = video ;
    end
end
%{
b = video;
b=reshape(b,size(video,1),size(video,2),size(video,3));
graysnap= rgb2gray(video) ;               %% converting color omage to grayscale

image(b);
image(video);
imshow(video);
imshow(graysnap);
% draw last pic
figure(size(snap,1));
%}


%a = reshape(snap(size(snap,1),:,:,:),size(snap,2),size(snap,3),size(snap,4));              % one  2D  frame
image(a);

b = rgb2gray(a);

imshow(b);
figure(2);

imshow(b);

level = 0.5 ;
Ithres = im2bw(b,level) ;       %% binary  segmantation  by threshold 'level'
Ithres = imbinarize(b);
imshow(Ithres,[]);
%{
% all image proccessing functions:
% https://www.mathworks.com/help/images/referencelist.html?type=function

Igray = rgb2gray(I) ;               %% converting color omage to grayscale
% probably better to take only the green image

level = 0.5 ;
Ithres = im2bw(Igray,level) ;       %% binary  segmantation  by threshold 'level'

J = imopen(I,SE)        ;           %% Morphologically open image
J = imclose(I,SE)       ;           %% Morphologically close image
%}