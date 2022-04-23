%% simulated_image_7b.m

%% creates a graphical simulation of one simulation results



folder = 'D:\Ben\EVOLVER\chromaffin endocytosis' ;


%path1 = 'E:\Ben\neuro-endocytosis\aug 2019 images';   %%output location
%{
movies = load(fullfile(folder,'movies')) ;
movie  = movies.movie;
str = strcat('choose movie number (1-',num2str(length(movie)),')');
video_number       = inputdlg(str); %% int8(str2num(inputdlg(str)));

scale =  movie(str2num(video_number{1})).scale_meter ;

snapshots_number = length(movie(str2num(video_number{1})).snapshots) ;
stru  = struct2cell(movie(str2num(video_number{1})).snapshots) ;

i = [7 8 9 ] ;  %% columns of h, ra and r0
I = [];
for counter = 1 : length(stru(7,1,:)) %% only where simulation location exists
   if not(isempty(stru{10,1,counter}))
    I = [I , counter];
   end
end

P  = cell2mat(stru(i,1,I));
P  = reshape(P,size(P,1),size(P,3));

r0  = P(1,:);
ra  = P(2,:);
H   = P(3,:);


Pn   = [H./r0 ; ra./r0] ;
%scatter(Pn(1,:),Pn(2,:));

kappa  = 0.8*10^(-19);
PA1   = r0.^2.*interp2(hq,rq,nAq,H./r0,ra./r0) ;    %% m^2 (interpolation is unitless)
%PA2   = PA1 - pi*ra.^2;                             %% m^2  -> the new
%version is after the subtraction
Pforce = (kappa./r0).*interp2(hq,rq,nforce,H./r0,ra./r0) ;                    %% points normal force [newton] = [kappa/r0]
Plambda= (kappa./r0).*interp2(hq,rq,nlambda,H./r0,ra./r0) ;                   %% points normal line tension [newton] = [kappa/r0]

%frames_index = dsearchn(E(:,[1,2],Pn) ;         % find the indices of simulation nearest to the fitted parametres.


% create the video writer with 1 fps
writerObj = VideoWriter('Simulated Video.avi');
writerObj.FrameRate = 1;
% set the seconds per image

open(writerObj); %% open the video for editing
%}
%a = load(movie(str2num(video_number{1})).snapshots(I(end)).simulation_location); %% look at final pic to determine axes limits
%a = a*movie(str2num(video_number{1})).snapshots(I(end)).natural_length_scale ;
%ax_limits = [-max(a(:,2)) max(a(:,2)) -max(a(:,2)) max(a(:,2)) min(a(:,1)) max(a(:,1))] ;

[files,path] = uigetfile('*.txt','Select One or More Files', 'MultiSelect', 'on');

max_h =0;
max_ra=0;
for snap_number = 1: length(files)  %% for finding the absolute largest dimension for setting fixed axes
    a = [];
    a = load(fullfile(path,files{snap_number}));
    %{
    %%%%%%%%%% if constant area is assumed !!
    %%%%%%%%%% scaling procedure %%%%%%%%%%%%
    AAA = interp2(hq,rq,tmpnAq, max(a(:,1)),a(2,2)) ;
    a = a*(Afixed/AAA).^0.5 ;
    %%%%%%%%%% scaling procedure %%%%%%%%%%%%
    %}
    if max_h < max(a(:,1))  %% for axes limits
        max_h = max(a(:,1));
    end
    if max_ra <  max(a(:,2))
        max_ra = max(a(:,2));
    end
end
ax_limits = [-max_ra max_ra -max_ra max_ra 0 max_h] ;
% loop that runs over all frames
% for each frame -> creats graphics, add description (i.e force) and saves an image
for snap_number = 1:length(files)
    
    %%%% graphics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    a = [];
    a =  load(fullfile(path,files{snap_number}));
   %{
        %%%%%%%%%% if constant area is assumed !!
        %%%%%%%%%% scaling procedure %%%%%%%%%%%%
        AAA = interp2(hq,rq,tmpnAq, max(a(:,1)),a(2,2)) ;
        a = a*(Afixed/AAA).^0.5 ;
        %%%%%%%%%% scaling procedure %%%%%%%%%%%%
        a(1,:) = [0 , max_ra];                                    %%get rid of the first vertex
    %}
    a(1,:) = [];  %%get rid of the first vertex
    % n= length(a(:,1));
    %{
        for i = 2:n                       %% bubble sort - arange matrix by x
            for j = 2:n-1 %$$$ weird (n-i)
                if a(j,1) > a(j+1,1)
                    tmp       = a(j,:);
                    a(j,:)    = a(j+1,:);
                    a(j+1,:)  = tmp;
                end
            end
        end
    %}
    
    [~,I1] = sort(a(2:end,1));                                                           % sort - arange matrix by z (skip first two)
    
    ZO  = a([1, I1(:)'+1],1);
    R   = a([1, I1(:)'+1],2);
    n   = length(ZO);
    
    resolution = 60;
    
    rmax = max(R);
    %rmax = max(a(:,2));
    %rmax = 10;
    
    
    dz =  min(abs(ZO-circshift(ZO,1)));
    if (dz==0)
        dz = 0.001;
    end
    
    m = floor((ZO(n)-ZO(1))/dz);
    
    %%% lets try->
    m = n;
    dz = (ZO(n)-ZO(1))/m;
    
    z0 = ZO(1);
    z=z0:dz:z0+(m)*dz;  % new z now equally spaced (for cylinder function)
    
    
    %{
        r=z;                            % tube's radius
        zleft = a(j,1);
        
        zright = a(j+1,1);
        for i=1:m                   %% make new equally distant z with interpolated r while (z(i)>zright)
        j=j+1;
        zleft = a(j,1);
        zright = a(j+1,1);
            end
           rleft = a(j,2);
           rright = a(j+1,2);
          r(i) = (rleft*(zright-z(i))+rright*(z(i)-zleft))/(zright-zleft);  %%simple linear interpolation
        end
    %}
    [zv,index] = unique(ZO(:));
    rv = R(index);
    r = interp1(zv,rv,z,'spline') ;  %% choose method
    % r = interp1(a(1:end-1,1),a(1:end-1,2),z,'spline');
    
    % r = [r,a(end,2)];
    % z = [z,a(end,1)];
    
    dilute = ceil(m/resolution);
    originalz = z;      %keep high resolution
    r = r(1:dilute:m);
    z = z(1:dilute:m);
    
    r = [r,R(end)];
    z = [z,ZO(end)];
    
    m = length(z(1,:));
    
    % r = a(:,2);
    [X,Y,Z] = cylinder(r);                  %% function the makes the meshgrid for surface --> r are rhe radiuses, z spreads in a unit
    
    % Z=z0+Z*(a(n,1)-a(1,1));        %% rescale z from unit to original dimensions
    Z=Z*(max(a(:,1))-min(a(:,1)));           %% rescale z from unit to original dimensions
    
    
    
    % color = 0.8*[255 225221 206]./255 ;  %% more colors: bez' - [0.7 0.6 0.5255] ; % redish - [0.7 0.46 0.46]
    color = [4 225 172]./255 ;  %% more colors: bez' - [0.7 0.6 0.5255] ; % redish - [0.7 0.46 0.46]
    figure('position',[500 150 600 800]);                               %position and most importently size of figure --> [500 150 600 800]
    surf(X,Y,Z,'faceColor',color,'Ambient',1.0);  %                                                                 $$$ opt3
    camlight('right')
    material dull;
    %surf(X,Y,Z);  %                                                                 $$$ opt3
    
    
    %surf(X,Y,Z,C,'Ambient',1.0);  %                                                                 $$$ opt3
    
    %h = light;
    %lighting gouraud;
    %axis  equal;
    %colormap pink;%parula;%winter;
    %caxis([0 100]);                                                         %%sets the bar limits of the specific colormap
    %axis([-1.5*rmax 1.5*rmax -1.5*rmax 1.5*rmax z0 -z0]);
    
    %  axis off;
    
    %grid off;
    
    %{
%% ring
hold on;
%[X,Y,Z] = cylinder([r(1) r(1)]);
[X,Y,Z] = cylinder(r(1));
Z=-0.1+Z*0.1;
surf(X,Y,Z,'Ambient',1.0);  %



%[X,Y,Z] = cylinder([r(end) 0.000]);
%Z= (Z./Z).*(z(end)+0.01); %z(end)+Z*0.01;
%surf(X,Y,Z,'Ambient',1.0);  %
viscircles([0 0 ],r(end));
circles(0,0,r(end),'color','black');
    %}
    
    
    %{
        %% plug
        hold on;
        color = [112, 69 46]./255;
        %[X,Y,Z] = cylinder([r(1) r(1)]);
        [X,Y,Z] = cylinder(r(end));
        %Z=z(end)+Z*0.01;
        Z=a(end,1)+Z*0.02;
        surf(X,Y,Z,'facecolor',color,'Ambient',1.0);  %
        %fill3(X(1,:),Y(1,:),Z(1,:),'k') ;
        fill3(X(2,:),Y(2,:),Z(2,:),color);
        
        %ring
        viscircles([0 0],r(1),'color',color);
        
        % fiber
        fiberend = 3.0;
        [X,Y,Z] = cylinder(r(end),6);
        Z=a(end,1)+Z*(fiberend-z(end));
        surf(X,Y,Z,'facecolor',color);%,'Ambient',1.0);  %
        %fill3(X(1,:),Y(1,:),Z(1,:),'k') ;
        %fill3(X(2,:),Y(2,:),Z(2,:),'r');
        
        %pos = [-r(end) -r(end) 2*r(end) 2*r(end)];
        %shape = rectangle('Position',pos,'Curvature',[1 1],'faceColor',[0 0 0]);
        %translate(shape,[0 0 z(end)]);
    %}
    
    axis equal ;
    axis (ax_limits) ;
    grid off ;
    %axis normal
    
    
    %%{
  %  set(gca,'XTick',[], 'YTick', [],'ZTick', [])
    axis off;
    %}
    
    %%%%%%% add description %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     str1 = join([num2cell('\bullet f =',2) num2cell(num2str(Pforce(snap_number==I)),2) num2cell('N',2)],2 );
    %     str2 = join([num2cell('\bullet \lambda =',2) num2cell(num2str(Plambda(snap_number==I)),2) num2cell('N',2)],2 );
    %     str3 = join([num2cell('\bullet A =',2) num2cell(num2str(PA1(snap_number==I)),2) num2cell('m^2',2)],2 );
    %     str = [str1 ; str2; str3 ] ;
    
    %text(min(min(X)),max(max(Y)),max(max(Z)),str) ;
    %     text(ax_limits(1),ax_limits(4),ax_limits(6),str) ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %save('try.bmp');
    %      status = mkdir(fullfile(folder,subfolder{folderindex}), 'illustrations') ;   %% creates new folder for pics
    %save(fullfile(folder,subfolder{folderindex},'illustrations',strcat('ra=',num2str(r(1)),'_h=',num2str(a(end,1)),'.bmp')));
    %     saveas(gcf,fullfile(folder,subfolder{folderindex},'illustrations',strcat('ra=',num2str(r(1)),'_h=',num2str(a(end,1)),'.jpg')));                                     %% save figure
    saveas(gcf,fullfile(path,strcat('snapshot-',num2str(snap_number),'.png')));                                     %% save figure
    
    
    
    %F = getframe(gcf) ;
    %writeVideo(writerObj, F);
    
    
    
    close all;
    
end
% close(writerObj);
% mov = immovie(writerObj)
%%%%%%%%%%%%%%%%%% make movie %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

