% shape_analysis_7b.m
% the curve representing shape transition
% r_a/r_0 vs h/r_0

%cd 'E:\Ben\EVOLVER\chromaffin endocytosis\harbe'
results_folder = 'D:\Ben\EVOLVER\chromaffin endocytosis\results august (mixes ra and h)' ;
cd(results_folder) ;


folders = dir('*');  %% 
dirFlags = [folders.isdir];
folders = folders(dirFlags);

r_a = [];
h= [];
surface_area = [];
tubeflag=[];
omegaflag=[];
lambdaflag=[];
f=[];
lamb=[];
Anorm=[];
A1norm=[];
A2norm=[];
stablef=[];
stablel=[];
neck_r = [];
threshold = tan(8*(pi/180));%10^-2; 8 degrees deviation
angle = [];

for folderindex = 1:length(folders)  % maybe start with 3
    cd(fullfile(results_folder,folders(folderindex).name));
    txt_files = dir('aperture*height*.txt');
    for count = 1 : length(txt_files)
        %% scan the shape from z=0 and stops upon reaching omega shape
        %count = 300;
        filename    = txt_files(count).name;
        % if exist(fullfile(folders(folderindex).folder,folders(folderindex).name,filename), 'file') == 2          %% check for existence of the .txt, if not, skip
        a = load(fullfile(folders(folderindex).folder,folders(folderindex).name,filename));
        
        
        %   if (max(z)>min(min(hq)))&&(max(z)<max(max(hq)))&&(r(2)>min(min(rq)))&&(max(z)<max(max(rq)))
        
        [~,I] = sort(a(3:end,1));                                                           % sort - arange matrix by z (skip first two)
        
        z = a([1, 2, I(:)'+2],1);
        r = a([1, 2, I(:)'+2],2);
        
        
        if (max(z)>min(min(hq)))&&(max(z)<max(max(hq)))&&(r(2)>min(min(rq)))&&(r(2)<max(max(rq)))
            
            tmp  = (circshift(r,-1,1)-circshift(r,1,1))./(circshift(z,-1,1)-circshift(z,1,1)) ; % symmetric derivative
            rtag = tmp(2:end-1);  %% changed from end-2
            %area1= trapz(z(2:end-1),2*pi*r(2:end-1).*(1+rtag.^2).^0.5);
            area1= trapz(z(2:end-1),2*pi*r(2:end-1).*(1+rtag.^2).^0.5)+pi*r(end)^2;  %% add plug area
            tubeflag   = [tubeflag,(-threshold < max(rtag) && max(rtag) < threshold)] ;         % logical flag indicating tube shape
            omegaflag  = [omegaflag,(max(rtag) > threshold)] ;                                  % logical flag indicating omega shape
            lambdaflag = [lambdaflag, ~tubeflag(end) && ~omegaflag(end)] ;
            r_a = [r_a,r(2)];
            h   = [h,max(z)];
            surface_area = [surface_area,area1] ;  %area1
            neck_r = [neck_r , neck_size(r,z) ];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % force and line tension for each point
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %[E_h,E_r] = gradient(Eq,dr,dh);    %% important to keep the spacing in both directions equal!
            %E_r = E_r/(2*pi);
            
            % Eq = griddata(E(:,1),E(:,2),E(:,3),rq,hq);
            %f    = [f, griddata(rq,rq,force,r(2),max(z)) ] ;
            %lamb = [lamb, griddata(rq,rq,lambda,r(2),max(z)) ] ;
            
            %Vq = interp2(X,Y,V,Xq,Yq)
            %     if (max(z)>min(min(hq)))&&(max(z)<max(max(hq)))&&(r(2)>min(min(rq)))&&(max(z)<max(max(rq)))
            fq = interp2(hq,rq,force,max(z),r(2));
            f    = [f, fq ] ;
            lambq = interp2(hq,rq,lambda,max(z),r(2));
            lamb = [lamb, lambq ] ;
            Anormq= interp2(hq,rq,nAq,max(z),r(2));
            Anorm = [Anorm , Anormq] ;
            A1normq= interp2(hq,rq,nAq+pi*rq.^2,max(z),r(2));
            A1norm = [A1norm , A1normq] ;
            
            A2normq= interp2(hq,rq,nA2q,max(z),r(2));
            A2norm = [A2norm , A2normq] ;
                        
            
            stablef = [stablef, interp2(hq,rq,stability_f,max(z),r(2))>0];
            stablel = [stablel, interp2(hq,rq,stability_lam,max(z),r(2))>0];
            
            angle = [angle, atan(max(rtag))];
        end
        
        
        %  end
    end
    
end



neck   = griddata(h,r_a,neck_r,hq,rq);


tmpangle = angle ;
tmph=h;
tmpr_a = r_a ;


for i =  length(angle) : -1 : 1
    j=1;
    while (j <= length(outliers_parameters))
        if r_a(i)==outliers_parameters(j)&&h(i)==outliers_parameters(j)
            tmpangle(i) = [] ;
            tmph(i)     = [] ;
            r_a(i)      = [] ;
            j = length(outliers_parameters)+1 ;
        else
            j = j+1 ;
        end
    end
end

nangle = griddata(tmph,tmpr_a,tmpangle,hq,rq);

%{
mesh(hq,rq,neck);
xlabel('H');
ylabel('r_{b}');
zlabel('neck radius');

nangle(50,251) = nan;
nangle(35,141) =nan;

nnangle = fillmissing(nangle,'pchip');
mesh(hq,rq,nnangle+pi/2,'LineStyle','none','faceColor','flat');
view(0,90)
C = colorbar ;
C.Limits    = [0 pi];
%C.Ticks    = [-pi/2 -pi/4 0 pi/4  ] ;
C.Ticks     = [0 pi/4  pi/2 3*pi/4 pi] ;
C.TickLabels= { '0' , '\pi/4' , '\pi/2' , '\pi3/4' , '\pi' };
C.FontSize = 14 ;
%C.Label.String= '\psi' ;
hL = ylabel(C,'\psi','Fontsize', 16);     
set(hL,'Rotation',0);
xlabel('$$\frac{H}{r_{i}}$$','interpreter','latex','fontsize',16,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{r_{b}}{r_{i}}$$','interpreter','latex','fontsize',16,'Fontweight','bold');%,'inerpreter','latex');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
path = 'E:\Ben\EVOLVER\chromaffin endocytosis\figures' ;
saveas(gcf,fullfile('E:\Ben\EVOLVER\chromaffin endocytosis\figures','angle diagram (2).fig'));
saveas(gcf,fullfile('E:\Ben\EVOLVER\chromaffin endocytosis\figures','angle diagram (2).jpeg'));

%map = [ zeros(101,1) , [1:-0.01:0]' , zeros(101,1)] ;
%map = [ zeros(101,1) , [0:0.01:1]' , [0:0.01:1]'] ;
map = [[ ones(101,1) ,  ones(101,1) , [0:0.01:1]' ]  ;  [  [1:-0.01:0]' , ones(101,1) , [1:-0.01:0]'] ] ;
map = [[ ones(101,1) ,  [0.9:0.001:1]' , [0:0.01:1]' ]  ;  [  [1:-0.01:0]' , [1:-0.001:0.9]' , [1:-0.01:0]'] ] ;
map = [[ ones(101,1) ,  [0.9:0.001:1]' , [0:0.01:1]' ]  ;  [  [1:-0.01:0]' , [1:-0.01:0]'] , [1:-0.001:0.9]'  ] ;
map = [[ ones(1001,1) ,  [0.9:0.0001:1]' , [0:0.001:1]' ]  ;  [  [1:-0.001:0]' , [1:-0.001:0]'] , [1:-0.0001:0.9]'  ] ;
colormap(map);
colormap('default');
xlim([0.05 12]);
ylim([0.0 4]);

%outliers_index      = [];
%outliers_parameters = [];
%%{
%% filter noise ->
apertures = unique(nangle(:,2));     %% create array of al apertures
%% loop and fuck it ->
for inti = 1:length(apertures)                      %% going through each aperture to remove outliers
    iap          = find(nangle(:,2)==apertures(inti));   %% find indices in E of the specific aperture
    [tmpi, ii]   = sort(nangle(iap,1));                  %% sort the indices by increasing heights -> iap(ii)=new indices of specific aperture aranged by increasing heights
    indices = iap(ii);
    [nangle(indices,3), out_i] = hampel(nangle(indices,3));            %% filter outlieres in the energy - note- it treats the distribution of h as they were equal (they are not always) (replaces noise by median of six points neighborhood)
   % E(iap(ii),3) = hampel(E(iap(ii),3),6,11);            %% filter outlieres in the energy - note- it treats the distribution of h as they were equal (they are not always) (replaces noise by median of six points neighborhood)
    %% set for hampel filter 6 neighbors and 11-> number of allowed standard deviations
    
  %  outliers_index      = [outliers_index , indices(out_i)'] ;  %% keep indices that was out filtered out (exclude them from fitting)
   % outliers_parameters = [outliers_parameters, E(indices(out_i),[1 2])' ] ;
end

heights = unique(nangle(:,1));     %% create array of al heights
%% loop and fuck it ->
for inti = 1:length(heights)                      %% going through each height to remove outliers
    ih          = find(nangle(:,1)==heights(inti));    %% find indices in E of the specific height
    [tmpi, ii]   = sort(nangle(ih,2));                  %% sort the indices by increasing apertures -> ih(ii)=new indices of specific aperture aranged by increasing heights
    indices = ih(ii);
    [nangle(indices,3), out_i] = hampel(nangle(indices,3));            %% filter outlieres in the energy - note- it treats the distribution of h as they were equal (they are not always) (replaces noise by median of six points neighborhood)
  %  E(iap(ii),3) = hampel(E(iap(ii),3),10,3);            %% filter outlieres in the energy - note- it treats the distribution of h as they were equal (they are not always) (replaces noise by median of six points neighborhood)
    %E(iap(ii),3)
    
    %outliers_index = [outliers_index , indices(out_i)'] ;
end



%}

cal_Aq = griddata(r_a,h,surface_area,rq,hq);
cal_nAq = fillmissing(cal_Aq,'pchip');
%surf(hq,rq,cal_Aq-pi*rq.^2);
cal_nAq = hampel(cal_nAq);
cal_nAq = hampel(cal_nAq')';

% create matrix of ra vs color
color = (lambdaflag+2*tubeflag+3*omegaflag)./(lambdaflag+tubeflag+omegaflag) ;
%{
[stability_f,~] = gradient(force,dr,dh);    %% important to keep the spacing in both directions equal!
[~,stability_lam] = gradient(lambda,dr,dh);    %% important to keep the spacing in both directions equal!

stability_f = stability_f.*(stability_f>0)./stability_f;
stability_lam = stability_lam.*(stability_lam>0)./stability_lam;
%}

%%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
%%% shape diagram by the two length parameters (r_a,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M       = zeros(3,length(color)); %% row or column in shape diagram
M(1,:)  = h ;
M(2,:)  = r_a ;
M(3,:)  = color;

scatter(M(1,:),M(2,:),20,M(3,:),'filled','square');  %% (x,y,size,color,filled,shape)
xlabel('$$\frac{H}{r_{i}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{r_{b}}{r_{i}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');


%}
%{
title('invaginations shapes');
text(3.5,1.25,'{\omega}');
text(2.3,1.82,'funnel');
text(3.5,1.25,'{\omega}');
text(2.3,1.82,'funnel');
picname = strcat('shape diagram','.png');                     %% name of figure to save
saveas(gcf,fullfile(folder,picname));
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% shape diagram by the height and  aperture(H,r_a) %%% interpolated and nice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
%[x1,y1]= ginput;                 %% the line seperating omega and tethers
xq1 = min(x1):0.01:max(x1);
yq1 = interp1(x1,y1,xq1,'spline');
yq1 = [ 0 yq1 yq1(end)];
xq1 = [ 0 xq1 12 ];

%[x2,y2]= ginput;                 %% the line seperating lambda and tethers
xq2 = min(x2):0.01:max(x2);
yq2 = interp1(x2,y2,xq2,'spline');
yq2 = [ 0 yq2];
xq2 = [ 0 xq2];
curves = [x1 ; y1 ; x2 ; y2 ];
figure(2);
%axis equal;
%daspect([max(xq1) max(yq2) 1]);
patch([xq1 xq1(end)], [yq1 0],[241 238 249]/255);   % Omega
patch([xq1 xq1(end) fliplr(xq2)], [yq1 max(y2) fliplr(yq2)],[248 242 218]/255);   % tether
patch([xq2 0], [yq2 yq2(end)],[222 242 248]/255);   % lambda

axis equal;
%xlim([0-0.01 max(xq1)+0.01]);
xlim([0-0.01 12+0.01]);
ylim([0-0.01 max(yq2)]+0.01);

set(gcf,'Position',[100 100 500 250]) ; 
xlabel('$$\frac{H}{r_{i}}$$','interpreter','latex','fontsize',15,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{r_{b}}{r_{i}}$$','interpreter','latex','fontsize',15,'Fontweight','bold');%,'inerpreter','latex');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
%xtickweight('bold');
%ytickweight('bold');
set(gca, 'Layer','top')
title('H-r_{a} shape diagram');

%%smooths lines:
%{
opengl hardware
set(0,'DefaultLineLineSmoothing','on');
set(0,'DefaultPatchLineSmoothing','on');
opengl('OpenGLLineSmoothingBug',1);
%}


%%% adding shapes curves %%%
hold on;
q = load('E:\Ben\EVOLVER\chromaffin endocytosis\results august (mixes ra and h)\simulation base 7b (34)\aperture_3.500000_height_2.000000.txt'); %lambda shape
q = q/8;
[q(:,1), ind]   = sort(q(:,1));
q(:,2)          = q(ind,2);
q(:,1)=q(:,1)-max(q(:,1));       %% position the top at (0,0)
orig = [3.5 2 ] ;
plot([-fliplr(q(:,2)) q(:,2)]+orig(2), [fliplr(q(:,1)) q(:,1)]+orig(1),'color','k');

q = load('E:\Ben\EVOLVER\chromaffin endocytosis\results august (mixes ra and h)\simulation base 7b (46)\aperture_3.500000_height_8.050000.txt'); %lambda shape
q = q/8;
[q(:,1), ind]   = sort(q(:,1));
q(:,2)          = q(ind,2);
q(:,1)=q(:,1)-max(q(:,1));       %% position the top at (0,0)
orig = [3.5 9 ] ;
plot([-fliplr(q(:,2)) q(:,2)]+orig(2), [fliplr(q(:,1)) q(:,1)]+orig(1),'color','k');

q = load('E:\Ben\EVOLVER\chromaffin endocytosis\results august (mixes ra and h)\simulation base 7b (23)\aperture_0.650000_height_3.200000.txt'); %lambda shape
q = q/8;
[q(:,1), ind]   = sort(q(:,1));
q(:,2)          = q(ind,2);
q(:,1)=q(:,1)-max(q(:,1));       %% position the top at (0,0)
orig = [1 6 ] ;
plot([-fliplr(q(:,2)) q(:,2)]+orig(2), [fliplr(q(:,1)) q(:,1)]+orig(1),'color','k');
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% shape diagram by the force and  aperture(f,r_a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
M       = zeros(3,length(color)); %% row or column in shape diagram
M(1,:)  = f ;
M(2,:)  = r_a ;
M(3,:)  = color;

scatter(M(2,:),M(1,:),20,M(3,:),'filled','square');  %% (x,y,size,color,filled,shape)
figure(2);
I=logical(color==1).*stablef;
scatter(M(2,logical(color==1).*stablef==1),M(1,logical(color==1).*stablef==1),20,M(3,logical(color==1).*stablef==1),'filled','square');  %% (x,y,size,color,filled,shape)
hold on;
line([0 max(M(2,I==1))], [2*pi 2*pi], 'color', 'red', 'linestyle', '--') ;
xlabel('$$Aperture/\sqrt{\frac{\kappa_{m}}{2\gamma}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{f}{\sqrt{2\kappa_{m}\gamma}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
set(get(gca,'YLabel'),'Position',[-0.3, 10, 0]);                                 %%y label position

%%% draw diagram
%I=logical(logical(color==1).*stablef);
M1 = M(1,I==1);
M2 = M(2,I==1);
M3 = M(3,I==1);
[M1,I]  = sort(M1,'descend') ;  % sort by adescending force (find maximum)
M2      = M2(I);
M3      = M3(I);
[M2,iM,~]   = unique(M2,'stable');
M1          = M1(iM);
M3          = M3(iM);
scatter(M2,M1,20,M3,'filled','square');  %% (x,y,size,color,filled,shape)
poly1 = [M2 ; M1 ] ;
%scatter(M2,M1,20,M3);


figure(3);
I=logical(color==2).*stablef;
scatter(M(2,I==1),M(1,I==1),20,M(3,I==1),'filled','square');  %% (x,y,size,color,filled,shape)
hold on;
line([0 max(M(2,I==1))], [2*pi 2*pi], 'color', 'red', 'linestyle', '--') ;

figure(4);
I=logical(color==3).*stablef;
scatter(M(2,I==1),M(1,I==1),20,M(3,I==1),'filled','square');  %% (x,y,size,color,filled,shape)
hold on;
line([0 max(M(2,I==1))], [2*pi 2*pi], 'color', 'red', 'linestyle', '--') ;

%%% draw diagram
M1 = M(1,I==1);
M2 = M(2,I==1);
M3 = M(3,I==1);
[M1,I]  = sort(M1,'ascend') ;  % sort by ascending force (find maximum)
M2      = M2(I);
M3      = M3(I);
[M2,iM,~]   = unique(M2,'stable');
M1          = M1(iM);
M3          = M3(iM);
scatter(M2,M1,20,M3,'filled','square');  %% (x,y,size,color,filled,shape)
poly2 = [M2 ; M1 ] ;
poly2(:,11) = [] ;          % outlier

%% unstable top part of omega
I=logical(color==3);
M1 = M(1,I==1);
M2 = M(2,I==1);
M3 = M(3,I==1);
[M1,I]  = sort(M1,'descend') ;  % sort by ascending force (find maximum)
M2      = M2(I);
M3      = M3(I);
[M2,iM,~]   = unique(M2,'stable');
M1          = M1(iM);
M3          = M3(iM);
scatter(M2,M1,20,M3,'filled','square');  %% (x,y,size,color,filled,shape)
poly3 = [M2 ; M1 ] ;

%%%%%%%%%%%%%%%%
% create diagram
%%%%%%%%%%%%%%%%
% poly 1- maximal force for lambda shape; poly 2- minimal force for omeg shape; poly 3- maximal force for omega
% shape (should be the same, maybe discard);
%{
figure(5);
hold on;
xq = 0.2:0.01:4;
yq = spline(poly1(1,:),poly1(2,:),xq);
area(xq,yq,0,'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',.3);
%set(gca, 'color', [1. 1. 0.7]);  % background color
fill([xq max(xq)], [yq max(yq)], 'y','FaceAlpha',.3,'EdgeAlpha',.3);  % also colors background
line([min(xq) max(xq)], [2*pi 2*pi], 'color', 'red', 'linestyle', '--') ;
xq = 0.2:0.01:max(poly2(1,:));
yq = spline(poly2(1,:),poly2(2,:),xq);
%area(xq,yq,2*pi,'FaceColor','r','FaceAlpha',.3,'EdgeAlpha',.3);
%hold off;
xq = [xq, max(xq), min(xq)];
yq = [yq, 2*pi, 2*pi];
fill(xq, yq, 'r','FaceAlpha',.5,'EdgeAlpha',.5);
axis([min(xq) 4 0 max(poly1(2,:))]);

xlabel('$$Aperture/\sqrt{\frac{\kappa_{m}}{2\gamma}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{f}{\sqrt{2\kappa_{m}\gamma}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
set(get(gca,'YLabel'),'Position',[-0.2, 25, 0]);                                 %%y label position
text(-0.09,2*pi+0.2,'2\pi\rightarrow','color','r','fontsize',12); %,'fontweight','bold'
title('f-r_{a} shape diagram');
hold off;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% shape diagram by the line tension and  height(lambda,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M       = zeros(3,length(color)); %% row or column in shape diagram
M(1,:)  = lamb ;
M(2,:)  = h ;
M(3,:)  = color;

figure(20);
scatter(M(2,:),M(1,:),20,M(3,:),'filled','square');  %% (x,y,size,color,filled,shape)
I=logical(stablel);
scatter(M(2,I==1),M(1,I==1),20,M(3,I==1),'filled','square');  %% (x,y,size,color,filled,shape)
figure(21);
I=logical(color==1).*stablel;
scatter(M(2,I==1),M(1,I==1),20,M(3,I==1),'filled','square');  %% (x,y,size,color,filled,shape)
hold on;
%line([0 max(M(2,I==1))], [2*pi 2*pi], 'color', 'red', 'linestyle', '--') ;
xlabel('$$Height/\sqrt{\frac{\kappa_{m}}{2\gamma}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{\lambda}{\sqrt{2\kappa_{m}\gamma}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
set(get(gca,'YLabel'),'Position',[-0.4, 8, 0]);                                 %%y label position
title('f-r_{a} shape diagram');
figure(31);
I=logical(color==2).*stablel;
scatter(M(2,I==1),M(1,I==1),20,M(3,I==1),'filled','square');  %% (x,y,size,color,filled,shape)
hold on;
%line([0 max(M(2,I==1))], [2*pi 2*pi], 'color', 'red', 'linestyle', '--') ;
figure(41);
I=logical(color==3).*stablel;
scatter(M(2,I==1),M(1,I==1),20,M(3,I==1),'filled','square');  %% (x,y,size,color,filled,shape)
hold on;
%line([0 max(M(2,I==1))], [2*pi 2*pi], 'color', 'red', 'linestyle', '--') ;



%%% draw diagram
I=logical(logical(color==1).*stablel);
M1 = M(1,I==1);
M2 = M(2,I==1);
M3 = M(3,I==1);
[M2,I]  = sort(M2,'ascend') ;  % sort by adescending height (find maximum)
M1      = M1(I);
M3      = M3(I);
[M1,iM,~]   = unique(M1,'stable');
M2          = M2(iM);
M3          = M3(iM);
scatter(M2,M1,20,M3,'filled','square');  %% (x,y,size,color,filled,shape)
poly4 = [M2 ; M1 ] ;
%scatter(M2,M1,20,M3);

%Omega:

%{
I=logical(logical(color==3).*stablel);
M1 = M(1,I==1);
M2 = M(2,I==1);
M3 = M(3,I==1);
[M1,I]  = sort(M1,'ascend') ;  % sort by ascending line tension (find maximum)
M2      = M2(I);
M3      = M3(I);
[M2,iM,~]   = unique(M2,'stable');
M1          = M1(iM);
M3          = M3(iM);
scatter(M2,M1,20,M3,'filled','square');  %% (x,y,size,color,filled,shape)
M1 = M1(M2>=0.11);
M3 = M3(M2>=0.11);
M2 = M2(M2>=0.11);
poly5 = [M2 ; M1 ] ;
%poly2(:,11) = [] ;          % outlier
%}


%%%%%%%%%%%%%%%%
% create diagram
%%%%%%%%%%%%%%%%
% poly 4- maximal line tension for lambda shape; poly 5- minimal line tension for omeg shape;
figure(30);
I=logical(color==1).*stablel;
scatter(M(2,I==1),M(1,I==1),20,M(3,I==1),'filled','square');  %% (x,y,size,color,filled,shape)
[x,y]= ginput;
poly4 = [x y];

% redefine poly4
figure(5);
hold on;
xq = 0.2:0.01:4;
%k  = convhull([poly4(1,:) max(poly4(1,:))],[poly4(2,:) max(poly4(1,:))]);
yq = interp1(poly4(:,1),poly4(:,2),xq,'cubic');
%poly4 = [poly4 ; [max(x) 0]];
%poly4 = [poly4 ; [min(x) 0]];
area(xq,yq,0,'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',.3);
%set(gca, 'color', [1. 1. 0.7]);  % background color
%fill([xq max(xq)], [yq max(yq)], 'y','FaceAlpha',.3,'EdgeAlpha',.3);  % also colors background

%{
% redefine poly5
figure(31);
I=logical(color==3).*stablel;
scatter(M(2,I==1),M(1,I==1),20,M(3,I==1),'filled','square');  %% (x,y,size,color,filled,shape)
[x,y]= ginput;
poly5 = [x y];

k  = convhull([poly5(1,:) max(poly5(1,:))],[poly5(2,:) max(poly5(1,:))]);
k(end)=[];
i = find(k==length(poly5(1,:))+1);
k(i)=[];
xq = 0.2:0.01:max(poly5(1,:));
[poly5x, i] = unique(poly5(1,k));
poly5y = poly5(2,k);
poly5y = poly5y(i);
%yq = spline(poly5x,poly5y,xq);
yq = interp1(poly5x,poly5y,xq,'pchip');
%area(xq,yq,2*pi,'FaceColor','r','FaceAlpha',.3,'EdgeAlpha',.3);
%hold off;
xq = [xq, max(xq)]; %, min(xq)
yq = [yq, max(yq)]; %, max(yq)

xq = 0.2:0.01:4;
%k  = convhull([poly4(1,:) max(poly4(1,:))],[poly4(2,:) max(poly4(1,:))]);
%yq = interp1(poly5(:,1),poly5(:,2),xq,'cubic');
yq = interp1(x,y,xq,'cubic');
poly4 = [poly4 ; [max(x) max(y)]];
poly4 = [poly4 ; [min(x) max(y)]];

xq = [xq, max(xq), min(xq)];
yq = [yq, max(yq), max(yq)];

%k  = convhull(xq,yq);
%fill(xq(k), yq(k), 'r','FaceAlpha',.5,'EdgeAlpha',.5);

%}

fill([xq, max(xq) , min(xq)], [yq, max(yq), max(yq)], 'r','FaceAlpha',.5,'EdgeAlpha',.5);
axis([min(xq) max(xq) 0 max(yq)]);
xlabel('$$Height/\sqrt{\frac{\kappa_{m}}{2\gamma}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{\lambda}{\sqrt{2\kappa_{m}\gamma}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
set(get(gca,'YLabel'),'Position',[-0.16, 5, 0]);
title('\lambda-H shape diagram');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% new diagram (f-lambda) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% with clamped area A0   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% shape diagram by the force and  line tension(f,lambda)%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pforce = interp2(hq,rq,nforce,1.5,4) ;                    %% points normal force 
Plambda= interp2(hq,rq,nlambda,1.5,4) ;                   %% points normal line tension 
PA     = interp2(hq,rq,nAq,1.5,4) ;                       %% points normal area

Afixed = PA;


M       = zeros(3,length(color));           %% row or column in shape diagram
M(1,:)  = f.*(Afixed./Anorm).^0.5 ;         %% newly  normalized force
M(2,:)  = lamb.*(Afixed./Anorm).^0.5 ;      %% newly  normalized line tension
M(3,:)  = color;

scatter(M(2,:),M(1,:),20,M(3,:),'filled','square');  %% (x,y,size,color,filled,shape)
I = stablef.*stablel;
scatter(M(2,I==1),M(1,I==1),20,M(3,I==1),'filled','square');  %% (x,y,size,color,filled,shape)

%[x20,y20]= ginput;                 %% the line seperating omega and tether
yq20 = min(y20):0.01:2*pi;
xq20 = interp1(y20,x20,yq20,'spline');

%[x21,y21]= ginput;                 %% the line seperating lambda and tether
xq21 = min(x20):0.01:max(xq21);
yq21 = interp1(x21,y21,xq21,'spline');
yq21 = min(yq21):0.01:2*pi;
xq21 = interp1(y21,x21,yq21,'spline');

%yq21(end) = 2*pi;

figure(8);
patch([fliplr(xq20) xq20(1) xq20(end)], [fliplr(yq20) 0 0],[241 238 249]/255);                  % Omega
patch([xq20  fliplr(xq21)], [yq20  fliplr(yq21)],[248 242 218]/255);                            % tether
patch([xq21 xq21(1)], [yq21 yq21(end)],[222 242 248]/255);                                      % lambda
patch([xq21(1) xq20(end) xq20(end) xq21(1)], [2*pi 2*pi 7 7],[202 202 202]/255);                % tether region


%axis equal;
xlim([min(xq20) max(xq20)]);
ylim([0-0.01 7]+0.01);

xlabel('$${\lambda}/\sqrt{\frac{\kappa_{m}}{2\gamma_{0}}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{f}{\sqrt{2\kappa_{m}\gamma_{0}}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
title('force-line tension shape diagram');


%%% adding shapes curves %%%
hold on;
q = load('E:\Ben\EVOLVER\chromaffin endocytosis\results august (mixes ra and h)\simulation base 7b (34)\aperture_3.500000_height_2.000000.txt'); %lambda shape
q = q./[6 80];
[q(:,1), ind]   = sort(q(:,1));
q(:,2)          = q(ind,2);
q(:,1)=q(:,1)-max(q(:,1));       %% position the top at (0,0)
orig = [5 0.3] ;
plot([-fliplr(q(:,2)) q(:,2)]+orig(2), [fliplr(q(:,1)) q(:,1)]+orig(1),'color','k');

q = load('E:\Ben\EVOLVER\chromaffin endocytosis\results august (mixes ra and h)\simulation base 7b (46)\aperture_3.500000_height_8.050000.txt'); %lambda shape
q = q./[11 147];
[q(:,1), ind]   = sort(q(:,1));
q(:,2)          = q(ind,2);
q(:,1)=q(:,1)-max(q(:,1));       %% position the top at (0,0)
orig = [5.3 0.95] ;
plot([-fliplr(q(:,2)) q(:,2)]+orig(2), [fliplr(q(:,1)) q(:,1)]+orig(1),'color','k');

q = load('E:\Ben\EVOLVER\chromaffin endocytosis\results august (mixes ra and h)\simulation base 7b (23)\aperture_0.650000_height_3.200000.txt'); %lambda shape
q = q./[6 80];
[q(:,1), ind]   = sort(q(:,1));
q(:,2)          = q(ind,2);
q(:,1)=q(:,1)-max(q(:,1));       %% position the top at (0,0)
orig = [1.05 0.95] ;
plot([-fliplr(q(:,2)) q(:,2)]+orig(2), [fliplr(q(:,1)) q(:,1)]+orig(1),'color','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for convexity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mesh(hq,rq,stability_f.*stability_lam)
view([0 0 1])
%[x31,y31]= ginput;                 %% bottom line
%[x32,y32]= ginput;                 %% top line
yq31 = min(y31):0.01:max(y31);
yq32 = min(y32):0.01:max(y32);
xq31 = interp1(y31,x31,yq31,'spline');
xq32 = interp1(y32,x32,yq32,'spline');

%yq31 = interp1(x31,y31,xq31,'linear');
xq31 = interp1(y31,x31,yq31,'linear');
%yq32 = interp1(x32,y32,xq32,'linear');
xq32 = interp1(y32,x32,yq32,'linear');

xq33 =[xq31 fliplr(xq32)];
yq33 =[yq31 fliplr(yq32)];

figure(2)
patch(xq33,yq33,[0.9,0.9,0.9]);

%{
picname = strcat('shape diagram','.png');                     %% name of figure to save
saveas(gcf,fullfile(folder,picname));
%}
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% shape diagram by the force and  aperture(f,r_a)  for constant area!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
Afixed = interp2(hq,rq,tmpnAq,2,3) ;
M       = zeros(3,length(color)); %% row or column in shape diagram
M(1,:)  = f.*(Afixed./abs(Anorm)).^0.5 ;
M(2,:)  = r_a.*(Afixed./abs(Anorm)).^0.5 ;
M(3,:)  = color;

scatter(M(2,:),M(1,:),20,M(3,:),'filled','square');  %% (x,y,size,color,filled,shape)

figure(2);
I=logical(color==1).*stablef;
scatter(M(2,logical(color==1).*stablef==1),M(1,logical(color==1).*stablef==1),20,M(3,logical(color==1).*stablef==1),'filled','square');  %% (x,y,size,color,filled,shape)
hold on;
line([0 max(M(2,I==1))], [2*pi 2*pi], 'color', 'red', 'linestyle', '--') ;
xlabel('$$Aperture/\sqrt{\frac{\kappa_{m}}{2\gamma_{0}}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{f}{\sqrt{2\kappa_{m}\gamma_{0}}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
set(get(gca,'YLabel'),'Position',[-0.3, 10, 0]);                                 %%y label position

%%% draw diagram
%I=logical(logical(color==1).*stablef);
M1 = M(1,I==1);
M2 = M(2,I==1);
M3 = M(3,I==1);
[M1,I]  = sort(M1,'descend') ;  % sort by adescending force (find maximum)
M2      = M2(I);
M3      = M3(I);
[M2,iM,~]   = unique(M2,'stable');
M1          = M1(iM);
M3          = M3(iM);
scatter(M2,M1,20,M3,'filled','square');  %% (x,y,size,color,filled,shape)
poly1 = [M2 ; M1 ] ;
%scatter(M2,M1,20,M3);


figure(3);
I=logical(color==2).*stablef;
scatter(M(2,I==1),M(1,I==1),20,M(3,I==1),'filled','square');  %% (x,y,size,color,filled,shape)
hold on;
line([0 max(M(2,I==1))], [2*pi 2*pi], 'color', 'red', 'linestyle', '--') ;

figure(4);
I=logical(color==3).*stablef;
scatter(M(2,I==1),M(1,I==1),20,M(3,I==1),'filled','square');  %% (x,y,size,color,filled,shape)
hold on;
line([0 max(M(2,I==1))], [2*pi 2*pi], 'color', 'red', 'linestyle', '--') ;

%%% draw diagram
M1 = M(1,I==1);
M2 = M(2,I==1);
M3 = M(3,I==1);
[M1,I]  = sort(M1,'ascend') ;  % sort by ascending force (find maximum)
M2      = M2(I);
M3      = M3(I);
[M2,iM,~]   = unique(M2,'stable');
M1          = M1(iM);
M3          = M3(iM);
scatter(M2,M1,20,M3,'filled','square');  %% (x,y,size,color,filled,shape)
poly2 = [M2 ; M1 ] ;
poly2(:,11) = [] ;          % outlier

%% unstable top part of omega
I=logical(color==3);
M1 = M(1,I==1);
M2 = M(2,I==1);
M3 = M(3,I==1);
[M1,I]  = sort(M1,'descend') ;  % sort by ascending force (find maximum)
M2      = M2(I);
M3      = M3(I);
[M2,iM,~]   = unique(M2,'stable');
M1          = M1(iM);
M3          = M3(iM);
scatter(M2,M1,20,M3,'filled','square');  %% (x,y,size,color,filled,shape)
poly3 = [M2 ; M1 ] ;

%%%%%%%%%%%%%%%%
% create diagram
%%%%%%%%%%%%%%%%
% poly 1- maximal force for lambda shape; poly 2- minimal force for omeg shape; poly 3- maximal force for omega
% shape (should be the same, maybe discard);
%{
figure(5);
hold on;
xq = 0.2:0.01:4;
yq = spline(poly1(1,:),poly1(2,:),xq);
area(xq,yq,0,'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',.3);
%set(gca, 'color', [1. 1. 0.7]);  % background color
fill([xq max(xq)], [yq max(yq)], 'y','FaceAlpha',.3,'EdgeAlpha',.3);  % also colors background
line([min(xq) max(xq)], [2*pi 2*pi], 'color', 'red', 'linestyle', '--') ;
xq = 0.2:0.01:max(poly2(1,:));
yq = spline(poly2(1,:),poly2(2,:),xq);
%area(xq,yq,2*pi,'FaceColor','r','FaceAlpha',.3,'EdgeAlpha',.3);
%hold off;
xq = [xq, max(xq), min(xq)];
yq = [yq, 2*pi, 2*pi];
fill(xq, yq, 'r','FaceAlpha',.5,'EdgeAlpha',.5);
axis([min(xq) 4 0 max(poly1(2,:))]);

xlabel('$$Aperture/\sqrt{\frac{\kappa_{m}}{2\gamma}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{f}{\sqrt{2\kappa_{m}\gamma}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
set(get(gca,'YLabel'),'Position',[-0.2, 25, 0]);                                 %%y label position
text(-0.09,2*pi+0.2,'2\pi\rightarrow','color','r','fontsize',12); %,'fontweight','bold'
title('f-r_{a} shape diagram');
hold off;
%}
%}