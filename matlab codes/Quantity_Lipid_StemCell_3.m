function []=Quantity_Lipid_StemCell_3(name_file,TRH_P_lipid,TRH_M_lipid,TRH_P_cell,TRH_M_cell,TRH_I,median_filter)
%apply the median filter directly on raw data nad then create the mask by
%using TRH of M and Intensity.
close all;
R64 = dir('*.R64'); 
%-------------------------------------------------------------------------------------------------------------------------------------
Nharm=1 ;% harmonic number
TT= 1000/(Nharm*80);% period expressed in ns that correspont to the laser frequancy repetition
om=2*pi/TT;
g_freeNADH= 1/(1+(om*0.4)^2);
s_freeNADH= om*0.4/(1+(om*0.4)^2);
tauref = 0:0.05:1000;
gdata=1./(1+(tauref*om).^2);
sdata=(tauref*om)./(1+(tauref*om).^2);
taudat=3.5;
g_3ns= 1/(1+(om*3)^2);
s_3ns=om*3/(1+(om*3)^2);
g_4ns= 1/(1+(om*4)^2);
s_4ns=om*4/(1+(om*4)^2);

%%%------------------------------------------------------------------------import raw data-------------------------------------------
%------512x512 pixel
k=length(R64);
for i = 1:length(R64)
T{i}=Open_R64();
I{i}=T{i}{1:1}(:,1:512);
I{i}=flip(rot90(T{1,i}{1:1}(:,1:512),3),2); 
P{i}=flip(rot90(T{1,i}{1:1}(:,513:1024),3),2);
M{i}=flip(rot90(T{1,i}{1:1}(:,1025:1536),3),2);

%%-----------------------------------------------------------create median filter---------------------------------------------------
Threshold_Inten{i}=I{i}>TRH_I;
I1{i}=medfilt2(I{i},median_filter);I1{i}=medfilt2(I1{i},median_filter);I1{i}=I1{i}.*Threshold_Inten{i};
P1{i}=medfilt2(P{i},median_filter);P1{i}=medfilt2(P1{i},median_filter);P1{i}=P1{i}.*Threshold_Inten{i};
M1{i}=medfilt2(M{i},median_filter);M1{i}=medfilt2(M1{i},median_filter);M1{i}=M1{i}.*Threshold_Inten{i};
g{i}=M1{i}.*cos(P1{i}.*pi/180);
s{i}=M1{i}.*sin(P1{i}.*pi/180);
%%--------------------------------------------------------------segmentation lipid droplets--------------------------------------------
Threshold_lipid_phase{i}=P1{i}>TRH_P_lipid;
Threshold_lipid_modul{i}=M1{i}<TRH_M_lipid;
mask_median_lipid{i}=Threshold_lipid_modul{i}.*Threshold_Inten{i}.*Threshold_lipid_phase{i};
Median_lipid_P{i}=P1{i}.*mask_median_lipid{i};
Median_lipid_M{i}=M1{i}.*mask_median_lipid{i};
tauphi_lipid{i}=(1/om).*(tan(Median_lipid_P{i}.*pi/180));
taumodul_lipid{i}=(1/om).*sqrt(1./Median_lipid_M{i}.^2-1);taumodul_lipid{i}(taumodul_lipid{i}==Inf)=0;
g_lipid{i}=Median_lipid_M{i}.*cos(Median_lipid_P{i}.*pi/180);
s_lipid{i}=Median_lipid_M{i}.*sin(Median_lipid_P{i}.*pi/180);
mean_tauP_lipid{i}=mean(mean(tauphi_lipid{i}(tauphi_lipid{i}>0)));
mean_tauM_lipid{i}=mean(mean(taumodul_lipid{i}(taumodul_lipid{i}>0)));

%%--------------------------------------------------------------segmentation mitochondria--------------------------------------------


Threshold_mito_phase{i}=P1{i}<TRH_P_cell;
Threshold_mito_modul{i}=M1{i}>TRH_M_cell;
mask_median_mito{i}=Threshold_Inten{i}.*Threshold_mito_modul{i}.*Threshold_mito_phase{i};

Median_cell_P{i}=P1{i}.*mask_median_mito{i};
Median_cell_M{i}=M1{i}.*mask_median_mito{i};
tauphi_mito{i}=(1/om).*(tan(Median_cell_P{i}.*pi/180));
taumodul_mito{i}=(1/om).*sqrt(1./Median_cell_M{i}.^2-1);taumodul_mito{i}(taumodul_mito{i}==Inf)=0;

g_cell{i}=Median_cell_M{i}.*cos(Median_cell_P{i}.*pi/180);%g and s of cell
s_cell{i}=Median_cell_M{i}.*sin(Median_cell_P{i}.*pi/180);%g and s of cell
mean_tauP_cell{i}=mean(mean(tauphi_mito{i}(tauphi_mito{i}>0)));
mean_tauM_cell{i}=mean(mean(taumodul_mito{i}(taumodul_mito{i}>0)));
%---------------------------------------------------------------------statistic------------------------------------------------------

%number of pixel between lipid and mitochondria 
pixel_lipid{i}=sum(sum(tauphi_lipid{i}(tauphi_lipid{i}~=0)));
pixel_mito{i}=sum(sum(tauphi_mito{i}(tauphi_mito{i}~=0)));
r{i}=pixel_lipid{i}/(pixel_mito{i}+ pixel_lipid{i})*100;
%-----lipid
 mean_g_lipid{i}=mean(mean(g_lipid{i}(g_lipid{i}~=0))); std_g_lipid{i}=std(g_lipid{i}(g_lipid{i}~=0));
 mean_s_lipid{i}=mean(mean(s_lipid{i}(s_lipid{i}~=0))); std_s_lipid{i}=std(s_lipid{i}(s_lipid{i}~=0));
%----cell
 mean_g_cell{i}=mean(mean(g_cell{i}(g_cell{i}~=0)));std_g_cell{i}=std(g_cell{i}(g_cell{i}~=0));
 mean_s_cell{i}=mean(mean(s_cell{i}(s_cell{i}~=0)));std_s_cell{i}=std(s_cell{i}(s_cell{i}~=0));
%--- calculate factor of bound NADH of mitochondria
%--------g and s of free NADH_tauphase=0.4ns
g_freeNADH= 1/(1+(om*0.4)^2);
s_freeNADH= om*0.4/(1+(om*0.4)^2);
fB_NADH{i}=sqrt(((g_cell{i}-g_freeNADH)).^2+((s_cell{i}-s_freeNADH)).^2);
mean_fB_NADH{i}=mean(mean(fB_NADH{i}(fB_NADH{i}<0.9)));
fB_NADH{i}(fB_NADH{i}>0.9)=0;
% plot mean g_s of every ROI
figure(4000);title('Phasor plot','fontsize',16);
plot(gdata,sdata,'-k','LineWidth',0.5);hold on;xlabel('g','fontsize',16);ylabel('s','fontsize',16);
plot(g_freeNADH,s_freeNADH,'--ks','LineWidth',2,'MarkerFaceColor','k','MarkerSize',4);label1 = {'0.4 ns'};
text(g_freeNADH,s_freeNADH,label1,'FontSize', 12,'VerticalAlignment','bottom','HorizontalAlignment','right');label4 = {'3 ns'};
plot(g_3ns,s_3ns,'--ks','LineWidth',2,'MarkerFaceColor','k','MarkerSize',4);%hold on;
text(g_3ns,s_3ns,label4,'FontSize', 12,'VerticalAlignment','bottom','HorizontalAlignment','right');label5 = {'4 ns'};
plot(g_4ns,s_4ns,'--ks','LineWidth',2,'MarkerFaceColor','k','MarkerSize',4);%hold on;
text(g_4ns,s_4ns,label5,'FontSize', 12,'VerticalAlignment','bottom','HorizontalAlignment','right');
errorbar(mean_g_lipid{i},mean_s_lipid{i},std_g_lipid{i}/2,std_g_lipid{i}/2, std_s_lipid{i}/2, std_s_lipid{i}/2,'.','MarkerSize',30,'MarkerEdgeColor','blue'); hold on;
errorbar(mean_g_cell{i},mean_s_cell{i},std_g_cell{i}/2,std_g_cell{i}/2, std_s_cell{i}/2, std_s_cell{i}/2,'.','MarkerSize',30,'MarkerEdgeColor','black'); hold on;
plot( mean_g_lipid{i},mean_s_lipid{i},'.','MarkerSize',30,'Color','blue');hold on;
plot( mean_g_cell{i},mean_s_cell{i},'.','MarkerSize',30,'Color','black');
xlim([0 1]);ylim([0 0.5]);hold on;legend('lipid','cell');hold on;
saveas (gcf,[char(name_file),'_g_s.fig']);saveas (gcf,[char(name_file),'_g_s.jpg']); 
%plot removed pixels
figure(i+22*k);
mask_remove_pixel_P{i}=(TRH_P_lipid>P1{i}).*(TRH_P_cell<P1{i}).*Threshold_Inten{i};
Removed_pixel_P_1{i}=mask_remove_pixel_P{i}.*P1{i}.*(TRH_P_lipid>P1{i}).*(TRH_P_cell<P1{i});
Removed_pixel_M_1{i}=mask_remove_pixel_P{i}.*M1{i}.*(TRH_P_lipid>P1{i}).*(TRH_P_cell<P1{i});
g_removed_pixel_1{i}=Removed_pixel_M_1{i}.*cos(Removed_pixel_P_1{i}.*pi/180);
s_removed_pixel_1{i}=Removed_pixel_M_1{i}.*sin(Removed_pixel_P_1{i}.*pi/180);
plot(g_removed_pixel_1{i},s_removed_pixel_1{i},'.','MarkerSize',1,'Color','black');hold on;
mask_remove_pixel_M{i}=(TRH_M_lipid<M1{i}).*(TRH_M_cell>M1{i}).*Threshold_Inten{i};
Removed_pixel_P_2{i}=mask_remove_pixel_M{i}.*P1{i}.*(TRH_M_lipid<M1{i}).*(TRH_M_cell>M1{i});
Removed_pixel_M_2{i}=mask_remove_pixel_M{i}.*M1{i}.*(TRH_M_lipid<M1{i}).*(TRH_M_cell>M1{i});
g_removed_pixel_2{i}=Removed_pixel_M_2{i}.*cos(Removed_pixel_P_2{i}.*pi/180);
s_removed_pixel_2{i}=Removed_pixel_M_2{i}.*sin(Removed_pixel_P_2{i}.*pi/180);
plot(g_removed_pixel_2{i},s_removed_pixel_2{i},'.','MarkerSize',1,'Color','black');hold on;
mask_remove_pixel_3{i}=(TRH_M_cell<M1{i}).*(TRH_P_lipid<P1{i}).*Threshold_Inten{i};
Removed_pixel_P_3{i}=mask_remove_pixel_3{i}.*P1{i}.*(TRH_M_cell<M1{i}).*(TRH_P_lipid<P1{i});
Removed_pixel_M_3{i}=mask_remove_pixel_3{i}.*M1{i}.*(TRH_M_cell<M1{i}).*(TRH_P_lipid<P1{i});
g_removed_pixel_3{i}=Removed_pixel_M_3{i}.*cos(Removed_pixel_P_3{i}.*pi/180);
s_removed_pixel_3{i}=Removed_pixel_M_3{i}.*sin(Removed_pixel_P_3{i}.*pi/180);
plot(g_removed_pixel_3{i},s_removed_pixel_3{i},'.','MarkerSize',1,'Color','black');hold on;
 mask_remove_pixel_4{i}=(TRH_M_lipid>M1{i}).*(TRH_P_cell<P1{i}).*Threshold_Inten{i};
Removed_pixel_P_4{i}=mask_remove_pixel_4{i}.*P1{i}.*(TRH_M_lipid>M1{i}).*(TRH_P_cell<P1{i});
Removed_pixel_M_4{i}=mask_remove_pixel_4{i}.*M1{i}.*(TRH_M_lipid>M1{i}).*(TRH_P_cell<P1{i});
g_removed_pixel_4{i}=Removed_pixel_M_4{i}.*cos(Removed_pixel_P_4{i}.*pi/180);
s_removed_pixel_4{i}=Removed_pixel_M_4{i}.*sin(Removed_pixel_P_4{i}.*pi/180);
 plot(g_removed_pixel_4{i},s_removed_pixel_4{i},'.','MarkerSize',1,'Color','black');hold on;
 mask_remove_pixel_5{i}=(TRH_P_lipid>P1{i}).*(TRH_M_cell>M1{i}).*Threshold_Inten{i};
Removed_pixel_P_5{i}=mask_remove_pixel_5{i}.*P1{i}.*(TRH_P_lipid>P1{i}).*(TRH_M_cell>M1{i});
Removed_pixel_M_5{i}=mask_remove_pixel_5{i}.*M1{i}.*(TRH_P_lipid>P1{i}).*(TRH_M_cell>M1{i});
g_removed_pixel_5{i}=Removed_pixel_M_5{i}.*cos(Removed_pixel_P_5{i}.*pi/180);
s_removed_pixel_5{i}=Removed_pixel_M_5{i}.*sin(Removed_pixel_P_5{i}.*pi/180);
 plot(g_removed_pixel_5{i},s_removed_pixel_5{i},'.','MarkerSize',1,'Color','black');hold on;
 %-----plot selected pixels
%-----plot cirlce of threshold
th = 0:pi/50:pi/2;
xunit1 = TRH_M_lipid * cos(th);
yunit1 = TRH_M_lipid * sin(th);
h1 = plot(xunit1, yunit1);hold on;%legend('TRH_M_lipid');
xunit2 = TRH_M_cell * cos(th);
yunit2 = TRH_M_cell * sin(th);
h2 = plot(xunit2, yunit2);hold on;%legend('TRH_M_cell');
xlim([0 1]);ylim([0 0.5]);
plot(gdata,sdata,'-k','LineWidth',0.5);hold on;xlabel('g','fontsize',16);ylabel('s','fontsize',16);
plot(g_freeNADH,s_freeNADH,'--ks','LineWidth',2,'MarkerFaceColor','k','MarkerSize',4);label1 = {'0.4 ns'};
text(g_freeNADH,s_freeNADH,label1,'FontSize', 12,'VerticalAlignment','bottom','HorizontalAlignment','right');label4 = {'3 ns'};
plot(g_3ns,s_3ns,'--ks','LineWidth',2,'MarkerFaceColor','k','MarkerSize',4);
text(g_3ns,s_3ns,label4,'FontSize', 12,'VerticalAlignment','bottom','HorizontalAlignment','right');label5 = {'4 ns'};
plot(g_4ns,s_4ns,'--ks','LineWidth',2,'MarkerFaceColor','k','MarkerSize',4);%hold on;
text(g_4ns,s_4ns,label5,'FontSize', 12,'VerticalAlignment','bottom','HorizontalAlignment','right');
plot( g_lipid{i},s_lipid{i},'.','MarkerSize',1,'Color','red');hold on;
plot( g_cell{i},s_cell{i},'.','MarkerSize',1,'Color','blue');hold on;
title(strcat('phasor plot ROI',num2str(i)),'fontsize',16);hold on;
saveas (gcf,[char(name_file),'_phasor plot.fig']);saveas (gcf,[char(name_file),'_phasor plot.jpg']);

end
%plot histogram of tauphi and taumodul
g_all=[g{i}];s_all=[s{i}];
[a,b]=size(g_all);
g_all=reshape(g_all,[a*b,1]);s_all=reshape(s_all,[a*b,1]);


phi_lipid_combine=[Median_lipid_P{i}];modul_lipid_combine=[Median_lipid_M{i}];
phi_cell_combine=[Median_cell_P{i}];modul_cell_combine=[Median_cell_M{i}];
 tauphi_lipid_combine= [tauphi_lipid{i}];
 tauphi_mito_combine= [tauphi_mito{i}];
 taumodul_lipid_combine= [taumodul_lipid{i}];
 taumodul_mito_combine= [taumodul_mito{i}]; 
 %calculate mean value of number of pixel between lipid and mitochondria
Values_1= [r].';T_1= table(Values_1);
writetable(T_1,[char(name_file) '_Ratio_lipid_mitochondria.txt']);
%---
Values_5= [mean_g_lipid].';T_5= table(Values_5);
writetable(T_5,[char(name_file) '_g_lipid.txt']);
%---
Values_6= [mean_s_lipid].';T_6= table(Values_6);
writetable(T_6,[char(name_file) '_s_lipid.txt']);
%---
Values_7= [mean_g_cell].';T_7= table(Values_7);
writetable(T_7,[char(name_file) '_g_cell.txt']);
%---
Values_8= [mean_s_cell].';T_8= table(Values_8);
writetable(T_8,[char(name_file) '_s_cell.txt']);
%---
Values_9= [mean_fB_NADH].';T_9= table(Values_9);
writetable(T_9,[char(name_file) 'mean_fB_NADH.txt']);
%---
ValuesNames_10={'TRH_P_lipid';'TRH_P_cell';'TRH_I';'TRH_M_lipid';'TRH_M_cell'};Values_10=[TRH_P_lipid;TRH_P_cell;TRH_I;TRH_M_lipid;TRH_M_cell];T_10= table(ValuesNames_10,Values_10);
writetable(T_10,[char(name_file) '_Threshold.txt']);
%----
Values_11= [mean_tauP_lipid].';T_11= table(Values_11);
writetable(T_11,[char(name_file) 'mean_tauP_lipid.txt']);
%----
Values_12= [mean_tauM_lipid].';T_12= table(Values_12);
writetable(T_12,[char(name_file) 'mean_tauM_lipid.txt']);
%----
Values_13= [mean_tauP_cell].';T_13= table(Values_13);
writetable(T_13,[char(name_file) 'mean_tauP_cell.txt']);
%----
Values_14= [mean_tauM_cell].';T_14= table(Values_14);
writetable(T_14,[char(name_file) 'mean_tauM_cell.txt']);
figure(200)
plot(gdata,sdata,'-k','LineWidth',0.5);hold on;xlabel('g','fontsize',16);ylabel('s','fontsize',16);
th = 0:pi/50:pi/2;
xunit1 = TRH_M_lipid * cos(th);
yunit1 = TRH_M_lipid * sin(th);
h1 = plot(xunit1, yunit1);hold on;%legend('TRH_M_lipid');
xunit2 = TRH_M_cell * cos(th);
yunit2 = TRH_M_cell * sin(th);
h2 = plot(xunit2, yunit2);hold on;%legend('TRH_M_cell');
xlim([0 1]);ylim([0 0.5]);
scatter_kde(g_all,s_all, 'filled', 'MarkerSize', 2);colormap(jet);
saveas (gcf,[char(name_file),'_phasor plot2.fig']);saveas (gcf,[char(name_file),'_phasor plot2.jpg']);
end
%end