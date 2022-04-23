% auto_fit_const_tension_7b.m

%fun = @(scale)mypoints(scale,q,z,r) ;            %% a is the fitting parameters ; t and f are the data vectors. my fun includes the calculation of A and B analytically


%folder = 'E:\Ben\EVOLVER\chromaffin endocytosis\results august (mixes ra and h)' ;
results_folder = 'E:\Ben\EVOLVER\chromaffin endocytosis\results august (mixes ra and h)' ;
moviefolder = 'E:\Ben\EVOLVER\chromaffin endocytosis' ;%% structure folder
path1 = 'E:\Ben\neuro-endocytosis\aug 2019 movies\movie 2 fit\constant tension';                     %% output folder
%folder = 'D:\Ben\EVOLVER\chromaffin endocytosis\results 30_01 U + 02_06' ;
%folder = 'D:\Ben\EVOLVER\chromaffin endocytosis\results 30_01 U + 02_06' ;
%folder = 'D:\Ben\EVOLVER\chromaffin endocytosis\results 30_01 U' ;
%files = dir('D:\Ben\EVOLVER\chromaffin endocytosis\results 30_01 U');
%%{
files = dir(results_folder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
%%}

alt_movies = load(fullfile(moviefolder,'alt_movies')) ;
%movies = load('D:\Ben\EVOLVER\chromaffin endocytosis\results 30_01 U + 02_06\movies') ;
alt_movie  = alt_movies.alt_movie;
str = strcat('choose movie number (1-',num2str(length(movie)),')');
video_number       = inputdlg(str); %% int8(str2num(inputdlg(str)));

mean_r_i = mean([movie(str2num(video_number{1})).snapshots(:).natural_length_scale]) ;
scale = mean_r_i/alt_movie(str2num(video_number{1})).scale_meter ;

number_of_folders_to_scan = 6 ;      % how many folders (of best averaged rough fit) to finally scan
test_number = 27;                    % how much simulations to check in first run in every folder

snapsss = 1 : 1031 ;
for snap_number  = (snapsss)
    leastsquares = -1;
    
    r   =   abs(alt_movie(str2num(video_number{1})).snapshots(snap_number).r);
    z   =   alt_movie(str2num(video_number{1})).snapshots(snap_number).z     ;
    rec =   alt_movie(str2num(video_number{1})).snapshots(snap_number).frame ;       %% frame of cropped image
    O   =   alt_movie(str2num(video_number{1})).snapshots(snap_number).axes_origin;
    alpha   =   alt_movie(str2num(video_number{1})).snapshots(snap_number).base_angle;
    
    
    
    folder_resnorm = zeros(1,length(subFolders)-2);  %% add 2 to the index! (first two folder are not folders
    for  folderindex = 3 : length(subFolders)   %% run over foldrs of simulations results
        cd(fullfile(results_folder,subFolders(folderindex).name));
        txt_files = dir('aperture*height*.txt');
        %% keep only data files in this folrder!
        %folder_resnorm=0;  %% the average value for goodness of fit
        
        leastsquares = -1 ;
        for filenumber = 1: floor(length(txt_files)/test_number) : length(txt_files)
            %filename = strcat(folder,'\',subFolders(k).name,'\',num2str(filenumber),'.txt');
            filename    = txt_files(filenumber).name;
            if isfile(filename)                          %% consition -> the file exists and not an outlier
                q = load(filename);
                aperture    = q(2,2);
                height      = max(q(:,1));
                outlier_check = ismember([height,aperture] , outliers_parameters, 'row' ) ;  %%  checking if the file is of a simulation that was an energy outlier
                if (outlier_check)   %% if true (an outlier) jump to next text file
                    continue
                end
               % scale0 = max(z)/height;
                 fun = @(scale)mypoints(scale,q,z,r) ;            %% a is the fitting parameters ; t and f are the data vectors. my fun includes the calculation of A and B analytically
                [scale, resnorm] = lsqnonlin(fun,scale0);
                
               
               
                if (resnorm < leastsquares) || (leastsquares == -1)
                    leastsquares    = resnorm ;
                end
            end
        end
        folder_resnorm(folderindex-2) = folder_resnorm(folderindex-2)+leastsquares ;  %% adds the lsq of best fit in folder (after running roughly)
    end
    [~,folder_rank] = sort(folder_resnorm(:));  %% sort folder by goodness of fit (add 2 to the index of "folder_index)
    
    folder_to_scan = folder_rank(1:number_of_folders_to_scan)+2 ;
    
    leastsquares = -1 ;
    for  folderindex = folder_to_scan'  % replaced '3 : length(subFolders)'   %% run over foldrs of simulations results
        cd(fullfile(results_folder,subFolders(folderindex).name));
        txt_files = dir('aperture*height*.txt');
        %% keep only data files in this folrder!
        for filenumber = 1:length(txt_files)
            %filename = strcat(folder,'\',subFolders(k).name,'\',num2str(filenumber),'.txt');
            filename    = txt_files(filenumber).name;
            if isfile(filename)                          %% consition -> the file exists and not an outlier
                
                q = load(filename);
                aperture    = q(2,2);
                height      = max(q(:,1));
                outlier_check = ismember([height,aperture] , outliers_parameters, 'row' ) ;  %%  checking if the file is of a simulation that was an energy outlier
                if (outlier_check)   %% if true (an outlier) jump to next text file
                    continue
                end
                scale0 = max(z)/height;
                fun = @(scale)mypoints(scale,q,z,r) ;            %% a is the fitting parameters ; t and f are the data vectors. my fun includes the calculation of A and B analytically
                [scale, resnorm] = lsqnonlin(fun,scale0);
                
                if (resnorm < leastsquares) || (leastsquares == -1)
                    leastsquares    = resnorm ;
                    bestfit         = filename ;
                    bestfitfolder   = subFolders(folderindex).name ;
                    bestscale       = scale ;
                    bestheight      = height ;
                    bestaperture    = aperture ;
                end
            end
        end
        
    end
    
    %cd 'E:\Ben\EVOLVER\chromaffin endocytosis';
    q = bestscale*load(fullfile(results_folder,bestfitfolder,bestfit));
    
    ROT = [cos(alpha) -sin(alpha) ; sin(alpha)  cos(alpha)]  ;  %% rotation matrix
    REF= [-1 0 ; 0 1]  ;  %% y-axis reflexion matrix
    %%{
    %q = load(bestfit);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ploting the fitted curve over the *original* picture
    %path1 = 'E:\Ben\neuro-endocytosis\jul 2019 movies\movie8';
    image(alt_movie(str2num(video_number{1})).snapshots(snap_number).og_image);
    hold on;
    [zfit, tmpind] = sort(q(:,1));
    rfit = q(tmpind,2);
    %plot(rec(1)+O(1)+zfit,rec(2)+O(2)+rfit);
    %plot(rec(1)+O(1)+rfit,rec(2)+O(2)+zfit);
    %plot(rec(1)+O(1)+rfit,rec(2)+O(2)-zfit);
    %figure(1)
    %hold on
    XY_p = [rec(1)+O(1) rec(2)+O(2)]' + ROT*[rfit  -zfit]' ;
    XY_m = [rec(1)+O(1) rec(2)+O(2)]' +ROT*REF*[rfit  -zfit]' ;
    
    plot(XY_p(1,:),XY_p(2,:),'r','LineWidth',2,'LineStyle',':');
    plot(XY_m(1,:),XY_m(2,:),'r','LineWidth',2,'LineStyle',':');
    %plot(rec(1)+O(1)+rfit,rec(2)+O(2)-zfit,'r','LineWidth',2,'LineStyle',':');
    %plot(rec(1)+O(1)-rfit,rec(2)+O(2)-zfit,'r','LineWidth',2,'LineStyle',':');
    title(strcat('snapshot- ',num2str(snap_number)));
    axis equal ;
    
    saveas(gca, fullfile(path1, strcat('fit over image snapshopt ',num2str(snap_number))), 'png');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    figure(2);  %% plot the fit over sample
    %scatter(movie(str2num(video_number{1})).scale_meter*q(:,2),movie(str2num(video_number{1})).scale_meter*q(:,1));
    [q(:,1),ind] = sort(q(:,1));
    q(:,2)       = q(ind,2);
    Z = [q(:,1)' fliplr(q(:,1)') ];
    R = [-q(:,2)' fliplr(q(:,2)')];
    %plot(movie(str2num(video_number{1})).scale_meter*q(:,2),movie(str2num(video_number{1})).scale_meter*q(:,1),'k');
    plot(alt_movie(str2num(video_number{1})).scale_meter*R,alt_movie(str2num(video_number{1})).scale_meter*Z,'k','LineWidth',1.5);
    hold on;
    r   =   (alt_movie(str2num(video_number{1})).snapshots(snap_number).r);
    Y =alt_movie(str2num(video_number{1})).scale_meter*[ z ; r ]';
    scatter(Y(:,2),Y(:,1),'*');
    title(strcat('snapshot- ',num2str(snap_number)));
    xlabel(strcat('r [m]'));
    ylabel(strcat('z [m]'));
    axis equal ;
    saveas(gca, fullfile(path1, strcat('fit snapshopt ',num2str(snap_number))), 'png');
    %[Idx, D] = knnsearch(q,Y);
    %Q = q(Idx,:);
    %scatter(Q(:,1),Q(:,2));
    %}
    
    alt_movie(str2num(video_number{1})).snapshots(snap_number).natural_length_scale  = bestscale*alt_movie(str2num(video_number{1})).scale_meter ;   %% calculating natural length scale
    alt_movie(str2num(video_number{1})).snapshots(snap_number).aperture = bestaperture*alt_movie(str2num(video_number{1})).snapshots(snap_number).natural_length_scale ;
    alt_movie(str2num(video_number{1})).snapshots(snap_number).height   = bestheight*alt_movie(str2num(video_number{1})).snapshots(snap_number).natural_length_scale   ;
    alt_movie(str2num(video_number{1})).snapshots(snap_number).simulation_location = fullfile(results_folder,bestfitfolder,bestfit) ;
    
    
    %save('D:\Ben\EVOLVER\chromaffin endocytosis\results 30_01 U + 02_06\movies','movie') ;
    cd 'E:\Ben\EVOLVER\chromaffin endocytosis' ;
    save('alt_movies','alt_movie') ;
    close all;
end


%save('movies','movie') ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% choosing for each Image point the nearest neighbor
function F = mypoints(scale,q,z,r)
%Y =scale*[ z ; r ]';
Y =[ z ; r ]';
[~, D] = knnsearch(scale*q,Y);
%Q = q(Idx,:);
F =    D;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%