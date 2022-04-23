% acquire_video.m
% what does it do?

folder      = 'D:\Ben\EVOLVER\chromaffin endocytosis\video' ;
filename    = 'media3.mp4';
%filename    = '03avi .mov';

v = VideoReader(fullfile(folder,filename));
numofframes = 8;
totalframes = floor(v.Duration*v.FrameRate);

%frames = read(v, 1:totalframes/numofframes:end);
counter = 0;
counter2 = 0;
snap = [];
while hasFrame(v)
    counter=counter+1 ;
    video = readFrame(v);
    if rem(counter,floor((totalframes)/(numofframes-1))) == 1
        counter2 = counter2+1;
        %onecolor = video(:,:,1)+video(:,:,2)+video(:,:,3);
        onecolor = video(:,:,2);
        snap(counter2,:,:) = onecolor;
        %figure(counter);
        %image(video);
    end
end



% draw last pic
figure(size(snap,1));
a = reshape(snap(size(snap,1),:,:),size(snap,2),size(snap,3));              % one  2D  frame
image(a);

%% set scale
realscalebar    = inputdlg('scale bar size in meters');
realscalebar    = str2double(realscalebar{1});
uiwait(helpdlg('click scalebar edges'));
[xi,yi]         = ginput(2);                                 % mark the scale bar
scalebar        = abs(xi(1)-xi(2));

%% crop image
uiwait(helpdlg('choose minimal window'));
[J , rec] = imcrop ;
rg = zeros(size(J,1),size(J,2));
close all;
%pic = [rg , J./max(J) , rg ];
%pic = reshape(pic,size(J,1),size(J,2),3);
%image(pic);


for index = size(snap,1):-1:1
    figure(index);
    a = reshape(snap(index,:,:),size(snap,2),size(snap,3));              % one  2D  frame
    %image(a);
    J = imcrop(a,rec) ;
    %image(J);
    pic = [rg , J./max(J) , rg ];
    pic = reshape(pic,size(J,1),size(J,2),3);
    image(pic);
    
    
    %[xi,yi]         = ginput(2);
    %line(xi,yi,'color','red');
    uiwait(helpdlg('set aperture edges'));
    answer = 'No';
    while answer(1) == 'N'
        p1 = impoint ;
        p2 = impoint ;
        p1pos = p1.getPosition;
        p2pos = p2.getPosition;
        x  = [p1pos(1) ; p2pos(1)];
        y  = [p1pos(2) ; p2pos(2)];
        apertureline    = line(x,y,'color','red');
        normal  = [y(2)-y(1),-(x(2)-x(1))];
        normal  = normal/(sqrt(normal*normal'));
        length = sqrt(size(pic,1).^2+size(pic,1).^2);
        %midline = line(([ sum(x)/2,sum(y)/2 ; sum(x)/2,sum(y)/2 ] + [-normal*length ; normal*length] ));
        midline = line([ sum(x)/2 sum(x)/2 ] + [-normal(1)*length ; normal(1)*length] , [ sum(y)/2 sum(y)/2 ] + [-normal(2)*length ; normal(2)*length]);
        uiwait(helpdlg('set tip point on the mid line'));
        p3 = impoint;
        answer = questdlg('the line is in the middle?');
        delete(p1);
        delete(p2);
        delete(p3);
        delete(apertureline);
        delete(midline);
    end
    
    
    uiwait(helpdlg('draw curve'));
    [xi,yi]         = ginput;
    
    answer = questdlg('to use this picture?');
    
    
end





