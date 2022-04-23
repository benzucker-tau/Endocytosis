%% surfaces_analysis_7b.m
%% interpolating surface using spline and smoothening
%% 
% energy vs r_a/r_0 & h/r_0
% E = (r_a,h,E)

cd 'D:\Ben\EVOLVER\chromaffin endocytosis\results august (mixes ra and h)'

%pulling_folders = dir('simulation*');
energies        = dir('**/Energy vs height*');


E =[];
A =[];
A2=[];
%% extracting vectors of aperture, height and energy from tables
for index = 1:length(energies)
    T = readtable(fullfile(energies(index).folder,energies(index).name));
    h1  = table2array(T(:,1));
    ra1 = table2array(T(:,2));
    e   = table2array(T(:,3));
    a   = table2array(T(:,4));
    a2  = table2array(T(:,5));
    if    isa(e,'cell')
    e=str2double(e);
    end
    h1(isnan(e))  = []; %% removing solutions with divergent energy
    ra1(isnan(e)) = [];
    a(isnan(e))   = [];
    a2(isnan(e))  = [];
    e(isnan(e))   = [];
    
    
   

    % needs to set the order
    %filtered_e = hampel(e);     %% filter outlieres in the energy - note- it treats the distribution of h as they were equal (they are not always) (replaces noise by median of six points neighborhood)
    
    %filtered_e = filloutliers(e,'spline');     %% filter outlieres in the energy - 
    %smoothed_e = smoothdata(filtered_e,'movmedian');     %% 'sgolay' 'movmedian'
    %smoothed_e = smooth(filtered_e);     %% 'sgolay' 'movmedian'
    
    %E = [E ;  h1, ra1, filtered_e];      %% E = ( r_a , h  , E)
    E = [E ;  h1, ra1, e];      %% E = ( r_a , h  , E)
    A = [A ; a] ;                               %% surface area for each point
    A2 = [A2 ; a2] ;                               %% second interpetation for surface area for each point (the whole area inside the aperture)
   % E = [E ; zeros(length(h1),1)+energies(folderindex).aperture, h1, smoothed_e];      %% E = ( r_a , h  , E)
end

%% bitchasss point I can't get rid of -> (try finding the problem)
%% manually eliminating ->


outliers_index      = [];
outliers_parameters = [];
%%{
%% filter noise ->
apertures = unique(E(:,2));     %% create array of al apertures
%% loop and fuck it ->
for inti = 1:length(apertures)                      %% going through each aperture to remove outliers
    iap          = find(E(:,2)==apertures(inti));   %% find indices in E of the specific aperture
    [tmpi, ii]   = sort(E(iap,1));                  %% sort the indices by increasing heights -> iap(ii)=new indices of specific aperture aranged by increasing heights
    indices = iap(ii);
    [E(indices,3), out_i] = hampel(E(indices,3));            %% filter outlieres in the energy - note- it treats the distribution of h as they were equal (they are not always) (replaces noise by median of six points neighborhood)
   % E(iap(ii),3) = hampel(E(iap(ii),3),6,11);            %% filter outlieres in the energy - note- it treats the distribution of h as they were equal (they are not always) (replaces noise by median of six points neighborhood)
    %% set for hampel filter 6 neighbors and 11-> number of allowed standard deviations
    
    outliers_index      = [outliers_index , indices(out_i)'] ;  %% keep indices that was out filtered out (exclude them from fitting)
   % outliers_parameters = [outliers_parameters, E(indices(out_i),[1 2])' ] ;
end

heights = unique(E(:,1));     %% create array of al heights
%% loop and fuck it ->
for inti = 1:length(heights)                      %% going through each height to remove outliers
    ih          = find(E(:,1)==heights(inti));    %% find indices in E of the specific height
    [tmpi, ii]   = sort(E(ih,2));                  %% sort the indices by increasing apertures -> ih(ii)=new indices of specific aperture aranged by increasing heights
    indices = ih(ii);
    [E(indices,3), out_i] = hampel(E(indices,3));            %% filter outlieres in the energy - note- it treats the distribution of h as they were equal (they are not always) (replaces noise by median of six points neighborhood)
  %  E(iap(ii),3) = hampel(E(iap(ii),3),10,3);            %% filter outlieres in the energy - note- it treats the distribution of h as they were equal (they are not always) (replaces noise by median of six points neighborhood)
    %E(iap(ii),3)
    
    outliers_index = [outliers_index , indices(out_i)'] ;
end


%filtered_e = hampel(e);     %% filter outlieres in the energy - note- it treats the distribution of h as they were equal (they are not always) (replaces noise by median of six points neighborhood)
%}

%{
%% outlier detection ->
apertures = unique(E(:,2));     %% create array of al apertures
%% loop and fuck it ->
for inti = 1:length(apertures)                      %% going through each aperture to remove outliers
    iap          = find(E(:,2)==apertures(inti));   %% find indices in E of the specific aperture
    [tmpi, ii]   = sort(E(iap,1));                  %% sort the indices by increasing heights -> iap(ii)=new indices of specific aperture aranged by increasing heights
    indices = iap(ii);
 %[TF,~,~,C] = isoutlier(E(iap(ii),3));
    B = isoutlier(E(indices,3),'median');
    %E(iap(ii),3) = B;
    %TF = isoutlier(A)
    %B = rmoutliers(A)
    %outliers_index = [outliers_index , indices(out_i)'] ;
    outliers_index = [outliers_index , indices(B)'] ;
end

heights = unique(E(:,1));     %% create array of al heights
%% loop and fuck it ->
for inti = 1:length(heights)                      %% going through each height to remove outliers
    ih          = find(E(:,1)==heights(inti));    %% find indices in E of the specific height
    [tmpi, ii]   = sort(E(ih,2));                  %% sort the indices by increasing apertures -> ih(ii)=new indices of specific aperture aranged by increasing heights
    indices = ih(ii);
    B = isoutlier(E(indices,3),'median');
    %E(iap(ii),3) = B;
    outliers_index = [outliers_index , indices(B)'] ;
end
%}

% fill in where NaN is present
%F = fillmissing(A,'constant',v);
%F = fillmissing(A,'constant',v);

eog=E;    %% keeps a copy of the original data
%{
f = 0;
lambda = 0; %%% 1.22121465;
E(:,3) = eog - f*E(:,2)+lambda*E(:,1);         %% F = U-h*f+r_a*lambda
%}

%% outliers_index = indices of simulations (points) filtered out
outliers_parameters = E(outliers_index,[1 2]) ;
E(outliers_index,:) = [];  %% remove outliers -> hampel filter interpolate energies but I still remove it and creates better interpolation
A(outliers_index) = [];  %% remove outliers 
A2(outliers_index) = [];  %% remove outliers 

figure(1)
scatter3(E(:,1),E(:,2),E(:,3),1,E(:,3));
xlabel('$$\frac{H}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{r_{a}}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
zlabel('$$energy/{\kappa}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
%hold on;

figure(2);
% smooth function interpolation
dr = 0.05;
dh = 0.05;
[~,min_ind] = min(E(:,1).^2+E(:,2).^2);
[~,max_ind] = max(E(:,1).^2+E(:,2).^2);  %% not a very good solution -> switch to a closer boundary for far boundary
[hq,rq] = meshgrid(E(min_ind,1):dh:E(max_ind,1),E(min_ind,2):dr:E(max_ind,2));
[hq,rq] = meshgrid(E(min_ind,1):dh:E(max_ind,1),E(min_ind,2):dr:5);
%[hq,rq] = meshgrid(min(E(:,1)):dr:max(E(:,1)),min(E(:,2)):dh:max(E(:,2)));
%[rq,hq] = meshgrid(min(r):dr:max(r),min(h1):dh:max(h1));
%
Eq = griddata(E(:,1),E(:,2),E(:,3),hq,rq,'cubic');
surf(hq,rq,Eq,'LineStyle',':','EdgeColor',[0.4,0.4,0.4]);%,'LineWidth',1.001);
xlabel('$$\frac{H}{r_{i}}$$','interpreter','latex','fontsize',15,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{r_{b}}{r_{i}}$$','interpreter','latex','fontsize',15,'Fontweight','bold');%,'inerpreter','latex');
zlabel('$$F_{el}/{\kappa}$$','interpreter','latex','fontsize',15,'Fontweight','bold');%,'inerpreter','latex');
xlim([0 12]);
ylim([0.1 4]);
saveas(gcf,fullfile('E:\Ben\EVOLVER\chromaffin endocytosis\figures','energy landscape.fig'));
saveas(gcf,fullfile('E:\Ben\EVOLVER\chromaffin endocytosis\figures','energy landscape.png'));

figure(7);
Eq = griddata(E(:,1),E(:,2),E(:,3),hq,rq);
%Eq = smoothdata(Eq,2);
%Eq = smoothdata(Eq,1);

%Eq = interp2(E(:,1),E(:,2),E(:,3),rq,hq,'spline');
%plot3(rq,hq,Eq);
surf(hq,rq,Eq);
%hold on;
xlabel('$$\frac{H}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{r_{a}}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
zlabel('$$energy/{\kappa}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');




%% force and line tension surfaces
figure(3);
%[force,lambda] = gradient(Eq);    %% important to keep the spacing in both directions equal!
[force,lambda] = gradient(Eq,dh,dr);    %% important to keep the spacing in both directions equal!
lambda = -lambda/(2*pi);
nforce = fillmissing(force,'pchip');
nlambda = fillmissing(lambda,'pchip');
[stability_f,stability_h_ra]   = gradient(force,dh,dr);    %% !
[stability_ra_h,stability_lam] = gradient(-lambda,dh,dr);    %% !
stability_f = stability_f.*(stability_f>0)./stability_f;
stability_lam = stability_lam.*(stability_lam>0)./stability_lam;
gauss_c = stability_f.*stability_lam-stability_h_ra.*stability_ra_h;
mean_c  = (stability_f+stability_lam)/2 ;
concaved = (gauss_c>0).*(mean_c>0);
mesh(hq,rq,concaved);
%contour(rq,hq,Eq);
%hold on;
%quiver(min(E(:,1)):dr:max(E(:,1)),min(E(:,2)):dh:max(E(:,2)),lambda,force);
mesh(hq,rq,force);
%hold off;
xlabel('$$\frac{H}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{r_{a}}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
set(get(gca,'YLabel'),'Position',[-0.31, 3, 0]);                                 %%y label position
zlabel('$$energy/{\kappa}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');

figure(4);
mesh(hq,rq,-lambda);
hold off;
xlabel('$$\frac{H}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{r_{a}}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
set(get(gca,'YLabel'),'Position',[-0.31, 3, 0]);                                 %%y label position
zlabel('$$energy/{\kappa}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');


%smf = smoothdata(force,1,'movmedian');     %% 'sgolay' 'movmedian'
%surf(rq,hq,-smf);

smf = smooth(force,2);
smf = reshape(smf,size(force,1),size(force,2));
surf(hq,rq,smf);

%plot(hq(1,:),gradient(Eq(:,10),dh));
%plot(hq(1,:),force(:,10))
%plot(hq(1,:),smooth(force(:,5)))
%plot(hq(1,:),smooth(smooth(smooth(force(:,5)))))
plot(hq(1,:),force(10,:));

%figure(4);
%%% first step streamlines - directed with energy gradient
%% written for old data type!
%{
start_ra = 1.0:1.0:3;
start_h = ones(size(start_ra))+3;
h = streamline(rq,hq,-lambda*2*pi,-force,start_ra,start_h); %% for streamline is more accurate to use lambda*2pi
paths_ra =  h.XData;
paths_h  =  h.YData;
paths_E  =   griddata(E(:,1),E(:,2),E(:,3),paths_ra,paths_h);
line(paths_ra,paths_h,paths_E,'color','r','LineWidth',3);

ylabel('$$\frac{H}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
xlabel('$$\frac{r_{a}}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
set(get(gca,'YLabel'),'Position',[-0.31, 3, 0]);                                 %%y label position
zlabel('$$energy/{\kappa}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
%}


figure(5);
%% area plot
Aq = griddata(E(:,1),E(:,2),A,hq,rq,'cubic');
nAq = fillmissing(Aq,'pchip');
A2q = griddata(E(:,1),E(:,2),A2,hq,rq,'cubic');
nA2q = fillmissing(A2q,'pchip');
surf(hq,rq,Aq);
%surf(hq,rq,nA2q);
%surf(hq,rq,Aq-Aq(:,1));
xlabel('$$\frac{H}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{r_{a}}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
zlabel('$$Area/{r_{0}^{2}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
title('surface area');

figure(6);
%% equi-forces contours
contour(hq,rq,force,[1  3  5  6 7 7.5 8],'ShowText','on','color','b');
%contour(hq,rq,force,[8.4-1.16*4 : 1.16 : 8.4],'ShowText','on','color','b');
hold on;
%line([0 3],[2.5 2.5] ) ;
contour(hq,rq,lambda,[0.1 0.3 0.5 1 1.1 1.2 1.3 1.4 1.5 2 2.5 4 6 ],'ShowText','on','color','g');
%patch(xq33,yq33,[0.9,0.9,0.9],'FaceAlpha',.5);
legend('puling force','constricting forces','unstable');
xlim([0 4.5]);
ylim([0 4.1]);

xlabel('$$\frac{H}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{r_{a}}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
title('equi-forces curves');

figure(7);
%% equi-forces contours for conserved area and varying tension for r_0 determined by the current tension
%% in this presentation one could apreciate the self-normalized parametrs of shapes for draw the shape and see structure (lambda, omega etc.)
nAq = fillmissing(Aq,'pchip');
nA2q = fillmissing(A2q,'pchip');
rb0 = 2.5 ;
h0  = 1;
%tmpnAq = nAq-nAq(:,1);  %% scale to get rid of negetive area -> could do it better
tmpnAq = nA2q; %% $$$$$$$$$$
%tmpnAq = nA2q;%-nAq(:,1);  %% scale to get rid of negetive area -> could do it better
%tmpnAq = nAq-min(min(nAq));  %% scale to get rid of negetive area -> could do it better

%Afixed = interp2(hq,rq,tmpnAq,5.5,3) ;                       %% normalized area of the shape when constricting begins ( relaxation with no membrane flow)
%Afixed = interp2(hq,rq,tmpnAq,1.2,3.5) ;                       %% (force=5.3) normalized area of the shape when constricting begins ( relaxation with no membrane flow)
%Afixed = interp2(hq,rq,tmpnAq,2,3) ;                       %% (force=5.3) normalized area of the shape when constricting begins ( relaxation with no membrane flow) ;                       %% (force=5.3) normalized area of the shape when constricting begins ( relaxation with no membrane flow)

Afixed = interp2(hq,rq,tmpnAq,h0,rb0) ;                       %% (force= 7.7) 
%Afixed = interp2(hq,rq,tmpnAq,1,2.5) ;                       %% (force=5.9) 

%Afixed = 2;
%normalization_c = (Afixed./nAq).^0.5;
%normalization   = real(normalization_c);
contour(hq,rq,force.*(Afixed./abs(tmpnAq)).^(-0.5),[4 6 interp2(hq,rq,nforce,h0,rb0) 8],'ShowText','on','color','b'); %% renormalized force (area conserved, tension changes)
hold on; 
%contour(hq,rq,lambda.*(Afixed./tmpnAq).^(-0.5),[interp2(hq,rq,lambda,1,2.5):(1.49-interp2(hq,rq,lambda,h0,rb0))/4  : 1.49],'ShowText','on','color','g');legend('puling force');
%xlim([0 5]);
%ylim([0 5]);
% contour(hq,rq.*(Afixed./abs(tmpnAq)).^0.5,force.*(Afixed./abs(tmpnAq)).^0.5,[1  3  5 5.3 5.6 5.9 7 7.5 8],'ShowText','on','color','b'); %% renormalized force (area conserved, tension changes)

contour(hq,rq,tmpnAq,[Afixed Afixed],'color','r'); %% renormalized force (area conserved, tension changes)

xlabel('$$\frac{H}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{r_{a}}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
title(strcat('equi-forces curves for A=',num2str(Afixed)));
% finding parameters from evenly spaced(self-choise) calculation-space parameters
%linelength = 1.5 ;
%line(linelength*[0 4.57],linelength*[0 0.5]) ;
%interp2(hq,nforce,tmpnAq,mat(:,1),mat(:,2))  

mat     = [1 2.5 ; 2.642 2.027 ; 3.912 1.51 ; 4.438 1.041 ; 4.553 0.5124]  ;                         %% real size
%mat     = [ 2.5 4.0000  ; 2.354 3.1000  ; 2.287  2.2000 ; 2.365 1.3000  ; 2.157  0.4000]  ;                         %% real size 4:-0.9:4-4*0.9
radius = 5 ;
angle  = atan(mat(:,2)./mat(:,1));
line([zeros(length(mat(:,1)),1) radius*cos(angle) ]',[zeros(length(mat(:,1)),1) radius*sin(angle) ]') ;  %% draw lines on normalized space to find shapes
normed = [1 2.5 ; 2.62 2.01 ; 4.82 1.86 ; 5.42 1.26 ; 5.42 0.61 ] ;
%normed = [ 2.5 4 ; 2.02 2.66 ; 1.77 1.71 ; 2.62 1.41 ; 2.92 0.56 ] ;
% normed = [1 2.5 ; 0.77 1.21 ; 5.22 1.56 ; 5.47 0.86 ; 5.37 0.51 ] ;     % for apcaced tensions
% for inserts:
% pore starts (two option)->  
% tether starts [3.02 2.11] (self-normed)
% force curve changes direction: [3.77 2.16]
% minimal rb  -> [5 0.15]   (self-normed)

%% ->   mat  = repmat((Afixed./interp2(hq,rq,tmpnAq,normed(:,1),normed(:,2))).^(0.5),1,2).*normed   %% self tension normalized
%forces check:
interp2(hq,rq,nforce,normed(:,1),normed(:,2)).*(Afixed./interp2(hq,rq,tmpnAq,normed(:,1),normed(:,2))).^(-0.5) 
interp2(hq,rq,nlambda,normed(:,1),normed(:,2)).*(Afixed./interp2(hq,rq,tmpnAq,normed(:,1),normed(:,2))).^(-0.5) 

figure(8);
%% equi-forces contours for conserved area and varying tension for the original r_0
%% in this presentation one could apreciate the true parametrs of shapes
rb0 = 2.5 ;
h0  = 1;
nAq = fillmissing(Aq,'pchip');
nA2q = fillmissing(A2q,'pchip');

%tmpnAq = nAq-nAq(:,1)+0.0001;  %% scale to get rid of negetive area -> could do it better

tmpnAq = nA2q;  %% scale to get rid of negetive area -> could do it better

%tmpnAq = nA2q;%-nAq(:,1);  %% scale to get rid of negetive area -> could do it better
%tmpnAq = nAq-min(min(nAq));  %% scale to get rid of negetive area -> could do it better

%Afixed = 2;
%Afixed = interp2(hq,rq,tmpnAq,2,3) ;                       %% (force=5.3) normalized area of the shape when constricting begins ( relaxation with no membrane flow)
Afixed = interp2(hq,rq,tmpnAq,h0,rb0) ;                       %% (force=5.3) normalized area of the shape when constricting begins ( relaxation with no membrane flow)
%normalization_c = (Afixed./nAq).^0.5;
%normalization   = real(normalization_c);
normalization   = (Afixed./abs(tmpnAq)).^0.5 ;
inv_normalization= normalization.^-1;
%contour(hq.*normalization,rq.*normalization,force.*inv_normalization,[   5 6 interp2(hq,rq,nforce,h0,rb0) 8],'ShowText','on','color','b'); %% renormalized force (area conserved, tension changes)
contour(hq.*normalization,rq.*normalization,force.*inv_normalization,[interp2(hq,rq,nforce,h0,rb0)  interp2(hq,rq,nforce,h0,rb0)],'ShowText','on','color','b'); %% renormalized force (area conserved, tension changes)
hold on;
%contour(hq.*normalization,rq.*normalization,lambda.*normalization,[0.5 1 1.1 1.2 1.3 1.4 1.5 2 2.5 4 6 ],'ShowText','on','color','g');
legend('puling force','constricting forces');
%xlim([0 5]);
%ylim([0 5]);

contour(hq.*normalization,rq.*normalization,tmpnAq,[Afixed Afixed],'color','r'); %% renormalized force (area conserved, tension changes)

xlim([0 6]);
ylim([0 4.5]);
xlabel('$$\frac{H}{r_{0(\gamma_{0})}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{r_{a}}{r_{0(\gamma_{0})}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
title(strcat('equi-forces curves for A=',num2str(Afixed)));
%% important note: the shape determines in the "self-tension" space and not in the "base-line" or "reservoir" normalized space

figure(9);
%% equi-area contours for conserved area and varying tension
%% nA2q or nAq
contour(hq,rq,nAq,[2:2:12],'ShowText','on','color','r'); %% renormalized force (area conserved, tension changes)
%xlim([0 5]);
%ylim([0 5]);

xlabel('$$\frac{H}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{r_{a}}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
title('equi-area curves');
%{
% save pic
picname = strcat('shape diagram','.png');                     %% name of figure to save
saveas(gcf,fullfile(folder,picname));

%}

