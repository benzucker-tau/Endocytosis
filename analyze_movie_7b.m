%% analyze_movie_7b.m

%% using the fitting results:
%% creates a path with time, h and ra.
%% draw the points on the parameters space
%% calculates the area of the shapes - A1 and A2 (minus the aperture disk)

folder = 'E:\Ben\EVOLVER\chromaffin endocytosis' ;
movies = load(fullfile(folder,'movies')) ; 
movie  = movies.movie;
str = strcat('choose movie number (1-',num2str(length(movie)),')');
video_number       = inputdlg(str); %% int8(str2num(inputdlg(str)));
%{
str = strcat('choose snapshot number (1-',num2str(size(movie(str2num(video_number{1})).snapshots,2)),')');
snap_number1       = inputdlg(str); %% int8(str2num(inputdlg(str)));
snap_number       = str2num(snap_number1{1}) ;% uint8(snap_number1{1}) ;
%}

scale =  movie(str2num(video_number{1})).scale_meter ;

snapshots_number = length(movie(str2num(video_number{1})).snapshots) ;
stru  = struct2cell(movie(str2num(video_number{1})).snapshots) ;

i = [8 9 10 ] ;  %% columns of h, ra and r0
I = [];
for counter = 1 : length(stru(7,1,:))
   if not(isempty(stru{7,1,counter}))
    I = [I , counter];   
   end
end
%I = I(1:115); % to cut movie 9
%{
qq = stru;
qqqq = reshape(qq,size(qq,1),size(qq,3));
qqq = qqqq(i,I) ;
qqqq = reshape(qqq,size(qqq,1),size(qqq,3));
%}
P  = cell2mat(stru(i,1,I));
P  = reshape(P,size(P,1),size(P,3));

r0  = P(1,:);
ra  = P(2,:);
H   = P(3,:);


Pn   = [H./r0 ; ra./r0] ;

scatter(Pn(1,:),Pn(2,:));

PA1   = r0.^2.*interp2(hq,rq,nAq,H./r0,ra./r0) ;    %% m^2 (interpolation is unitless)
PA2   = PA1 - pi*ra.^2;                             %% m^2
PA3   = r0.^2.*interp2(hq,rq,nA2q,H./r0,ra./r0) ;    %%

kappa  = 0.8*10^(-19);
Pforce = (kappa./r0).*interp2(hq,rq,nforce,H./r0,ra./r0) ;                    %% points normal force [newton] = [kappa/r0]
Plambda= (kappa./r0).*interp2(hq,rq,nlambda,H./r0,ra./r0) ;                   %% points normal line tension [newton] = [kappa/r0]
gamma  = 0.5*kappa./(r0.^2) ;
Pneck  = r0.*interp2(hq,rq,neck,H./r0,ra./r0) ;
%% scatter text at points position -> switch to area maybe?

%%%%% print area (PA1)
%text([H./r0],[ra./r0],(num2cell(PA1)));
figure(1);
str1 =  num2str(I');
str2 = repmat('\bullet [',length(I),1);
str3 = repmat('] A[m^{2}]=',length(I),1);
str = [str2 str1 str3] ;
%str2= strcat('{\n}'); %snapshots_number
str = join([num2cell(str,2) num2cell(num2str(PA1'),2) num2cell(num2str(PA1'),2)],2) ;
%str = join([num2cell(str1,2) num2cell(num2str(PA1'),2) num2cell(num2str(PA1'),2)],2) ;
%str = join([repmat({'\bullet'},4,1) num2cell(num2str(PA1'),2)],2)
text([H./r0],[ra./r0],str,'FontSize',8);
%xlim([min(H./r0) max(H./r0)]);
%ylim([min(ra./r0) max(ra./r0)]);
xlim([0 1.2*max(H./r0)]);
ylim([0 1.2*max(ra./r0)]);



figure(2);
str1 =  num2str(I');
str2 = repmat('\bullet [',length(I),1);
str3 = repmat('] A[m^{2}]=',length(I),1);
str = [str2 str1 str3] ;
str = join([num2cell(str,2) num2cell(num2str(PA2'),2)],2);
%str = join([repmat({'\bullet'},4,1) num2cell(num2str(PA1'),2)],2)
text([H./r0],[ra./r0],str,'FontSize',8);
%xlim([min(H./r0) max(H./r0)]);
%ylim([min(ra./r0) max(ra./r0)]);
xlim([0 1.2*max(H./r0)]);
ylim([0 1.2*max(ra./r0)]);
xlabel('$$\frac{H}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$\frac{r_{a}}{r_{0}}$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
set(get(gca,'YLabel'),'Position',[-0.31, 3, 0]);                                 %%y label position

%% print also pulling force and line tension
%{
figure(3);
str1= strcat('\bullet [', num2str(I),'] A[m^{2}]='); %(1:snapshots_number)'
str = join([num2cell(str1,2) num2cell(num2str(PA2'),2)],2)
%str = join([repmat({'\bullet'},4,1) num2cell(num2str(PA1'),2)],2)
text([H./r0],[ra./r0],str);
%xlim([min(H./r0) max(H./r0)]);
%ylim([min(ra./r0) max(ra./r0)]);
xlim([0 1.2*max(H./r0)]);
ylim([0 1.2*max(ra./r0)]);
%}

%% area plot over the real size (not normalized)
figure(4);
str1 =  num2str(I');
str2 = repmat('\bullet [',length(I),1);
str3 = repmat('] A[m^{2}]=',length(I),1);
str = [str2 str1 str3] ;
str = join([num2cell(str,2) num2cell(num2str(PA2'),2)],2);
%str = join([repmat({'\bullet'},4,1) num2cell(num2str(PA1'),2)],2)
text([H],[ra],str,'FontSize',8);
%xlim([min(H./r0) max(H./r0)]);
%ylim([min(ra./r0) max(ra./r0)]);
xlim([0 1.2*max(H)]);
ylim([0 1.2*max(ra)]);
xlabel('$${H}[m]$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$r_{a} [m]$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
set(get(gca,'YLabel'),'Position',10^-7*[-0.61, 2.8, 0]); 

%% force plot over the real size (real and not normalized)
figure(5);
str1 =  num2str(I');
str2 = repmat('\bullet [',length(I),1);
str3 = repmat('] f[N]=',length(I),1);
str = [str2 str1 str3] ;
str = join([num2cell(str,2) num2cell(num2str(Pforce'),2)],2);
%str = join([repmat({'\bullet'},4,1) num2cell(num2str(PA1'),2)],2)
text([H],[ra],str,'FontSize',8);
%xlim([min(H./r0) max(H./r0)]);
%ylim([min(ra./r0) max(ra./r0)]);
xlim([0 1.2*max(H)]);
ylim([0 1.2*max(ra)]);
xlabel('$${H}[m]$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$r_{a} [m]$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
set(get(gca,'YLabel'),'Position',10^-7*[-0.61, 2.8, 0]); 

%% line tension plot over the real size (real and not normalized)
figure(6);
str1 =  num2str(I');
str2 = repmat('\bullet [',length(I),1);
str3 = repmat('] {\lambda}[N]=',length(I),1);
str = [str2 str1 str3] ;
str = join([num2cell(str,2) num2cell(num2str(Plambda'),2)],2);
%str = join([repmat({'\bullet'},4,1) num2cell(num2str(PA1'),2)],2)
text([H],[ra],str,'FontSize',8);
%xlim([min(H./r0) max(H./r0)]);
%ylim([min(ra./r0) max(ra./r0)]);
xlim([0 1.2*max(H)]);
ylim([0 1.2*max(ra)]);
xlabel('$${H}[m]$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
ylabel('$$r_{a} [m]$$','interpreter','latex','fontsize',13,'Fontweight','bold');%,'inerpreter','latex');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
set(get(gca,'YLabel'),'Position',10^-7*[-0.61, 2.8, 0]); 

close all ;

%% plot time vs things
%%time vs force
figure (7);
%snaps = [ 6 44 141   251 375 ];   %% movie 1 snapshots
snaps = [ 6 44 102  251 375 ];   %% movie 1 snapshots
%snaps = [19  247 400 541  1025];   %% movie 2 snapshots
%snaps = [9 31 53 83 110];   %% movie 3 snapshots

scatter(I/movie(str2num(video_number{1})).frame_rate,Pforce*10^12,'s','filled','k','SizeData',10);
%scatter(I/movie(str2num(video_number{1})).frame_rate,Pforce*10^12,'d','k','SizeData',20);
%scatter(I/movie(str2num(video_number{1})).frame_rate,Pforce*10^12,'.','k','SizeData',20);
hold on;
scatter(I(snaps)/movie(str2num(video_number{1})).frame_rate,Pforce(snaps)*10^12,'s','filled','r','SizeData',20);
xlabel('$$t$$ [sec]','interpreter','latex','fontsize',15,'Fontweight','bold') ; 
ylabel('$$f$$$$[pN]$$','interpreter','latex','fontsize',15,'Fontweight','bold') ;
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
title('flat to \Lambda transition movie');
%title('\Lambda to \Omega transition movie');
%set(gcf,'Position',[100 100 280 210]) ;  % half - original pos = [~ ~ 560 420]
set(gcf,'Position',[100 100 300 210]) ;
ylim([0 5]);

figure (8);
scatter(I/movie(str2num(video_number{1})).frame_rate,Plambda*10^12);
xlabel('t [sec]');
ylabel('{\lambda} [pN]');
title('flat to \Lambda transition movie');
%title('\Lambda to \Omega transition movie');


figure (9);
scatter(I/movie(str2num(video_number{1})).frame_rate,ra*10^9,'s','filled','k','SizeData',10);
hold on;
scatter(I(snaps)/movie(str2num(video_number{1})).frame_rate,ra(snaps)*10^9,'s','filled','r','SizeData',20);
xlabel('$$t$$ [sec]','interpreter','latex','fontsize',15,'Fontweight','bold') ; 
ylabel('$$r_{b}$$$$[nm]$$','interpreter','latex','fontsize',15,'Fontweight','bold') ;
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
%title('flat to \Lambda transition movie');
%xlim([0 38]);
title('\Lambda to \Omega transition movie (3)');
set(gcf,'Position',[100 100 300 210]) ;
xlim([0 43]);


figure (9); % rb normalized
scatter(I/movie(str2num(video_number{1})).frame_rate,ra./r0,'s','filled','k','SizeData',20);
hold on;
scatter(I(snaps)/movie(str2num(video_number{1})).frame_rate,ra(snaps)./r0(snaps),'s','filled','r','SizeData',40);
xlabel('$$t$$ [sec]','interpreter','latex','fontsize',15,'Fontweight','bold') ; 
ylabel('$$\frac{r_{b}}{r_{i}}$$','interpreter','latex','fontsize',15,'Fontweight','bold') ;
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
%title('flat to \Lambda transition movie');
title('\Lambda to \Omega transition movie (2)');


figure (11);
scatter(I/movie(str2num(video_number{1})).frame_rate,gamma,'s','filled','k','SizeData',20);
hold on;
scatter(I(snaps)/movie(str2num(video_number{1})).frame_rate,gamma(snaps),'s','filled','r','SizeData',40);
xlabel('$$t$$ [sec]','interpreter','latex','fontsize',15,'Fontweight','bold') ; 
ylabel('{\gamma} [\muN/n]');
title('flat to \Lambda transition movie');
%title('\Lambda to \Omega transition movie');

scatter(I/movie(str2num(video_number{1})).frame_rate,r0*10^6,'s','filled','k','SizeData',20);
%hold on;
%scatter(I(snaps)/movie(str2num(video_number{1})).frame_rate,r0(snaps),'s','filled','r','SizeData',40);
xlabel('$$t$$ [sec]','interpreter','latex','fontsize',15,'Fontweight','bold') ; 
ylabel('$${r_{i}}$$ [$$\mu$$m]','interpreter','latex','fontsize',15,'Fontweight','bold') ;
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle');
%title('flat to \Lambda transition movie');
title('\Lambda to \Omega transition movie');
hold on ;
plot([0 20] , [1 1]*mean(r0(20*12:end))*10^6,'r--');
%title('\Lambda to \Omega transition movie');



figure (12);
scatter(I/movie(str2num(video_number{1})).frame_rate,PA2);
xlabel('t [sec]');
%ylabel('{\gamma} [\muN/n]');
title('flat to \Lambda transition movie');

figure (13);
scatter(I/movie(str2num(video_number{1})).frame_rate,Pneck);
xlabel('t [sec]');
xlim([0 max(I/movie(str2num(video_number{1})).frame_rate)]);
ylabel('neck [m]');
title('\Lambda to \Omega transition movie');
%title('flat to \Lambda transition movie');

figure(14) ; 
%% present sequence of shapes from movie on a time scale
%% show time developing series. top row green + dashed fit, bottom row scatter + countinuous fit

frames = [ 1 ];  % frames to show




axes('pos',[.1 .6 .5 .3]) ;   %% new axes for image -> [bottomleftcornerXposition bottomleftcornerYposition width height]
imshow('coins.png') ;



%scatter(I,Pforce);
%scatter(I,PA3);
%scatter(I,gamma);
%{
mean(fillmissing(Pforce,'pchip'))
std(fillmissing(Pforce,'pchip'))

mean(fillmissing(Plambda,'pchip'))
std(fillmissing(Plambda,'pchip'))

mean(fillmissing(gamma,'pchip'))
std(fillmissing(gamma,'pchip'))

mean(fillmissing(PA3,'pchip'))
std(fillmissing(PA3,'pchip'))
%}


