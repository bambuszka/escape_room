%% Code computes the p-shift motifs given a time-ordered edgelist of interactions 
%% such that the maximum gap of time (in seconds) allowed between successive interactions for them to be considered as a p-shift motif
THRESHOLD_TIME_GAP=5 % 

%% INPUT vaiables decription:
% [ts,tf, n1,n2 ,signs] : [start time of interaction, end time of interaction, source node id, target node id, emotional sign of interaction]


%%
count=0;close all
TeamIds=zeros(40,1);
for i=1:65
filename=strcat('team',...
num2str(i),'_updatedSubres_Complex.txt');
if isfile(filename)==1
    count=count+1;
    TeamIds(count)=i;
end
end

%% Storing Succesful and Failed Teams
countFail=0;
countSuccess=0;
TIMETAKEN=zeros(40,1);
indF=zeros(26,1);
indS=zeros(14,1);
for id=1:40
filename=strcat('team',...
        num2str(TeamIds(id)),'_updatedSubres_Complex.txt');
[ts,tf, n1,n2 ,signs]= readvars(filename);
ind=find(~isnan(ts));
tsf=ts(ind);tff=tf(ind);n1f=n1(ind);n2f=n2(ind);signsf=signs(ind);

durations=seconds(tf-ts);

Tf=max(tff);Ti=min(tsf);
TIMETAKEN(id)=hours(Tf-Ti);
if TIMETAKEN(id)>1
countFail=countFail+1;   
    indF(countFail)=id;
else
    countSuccess=countSuccess+1;
    indS(countSuccess)=id;
    
end
end

[TeamIds, TIMETAKEN];
%%
t=TIMETAKEN>1;TIMETAKEN60=TIMETAKEN.*(1-t)+t;

%% Removing data above 60 mins
%% and counting motifs
store_motif_times_S_null=zeros(13,1);count_motif_S_null=zeros(13,1);
store_motif_times_F_null=zeros(13,1);count_motif_F_null=zeros(13,1);
store_motif_times_S=zeros(13,1);count_motif_S=zeros(13,1);
store_motif_times_F=zeros(13,1);count_motif_F=zeros(13,1);


for id=1:40
    id
filename=strcat('team',...
        num2str(TeamIds(id)),'_updatedSubres_Complex.txt');
[ts,tf, n1,n2 ,signs]= readvars(filename);
ind1=find(~isnan(ts));ind2=find( minutes(tf(ind1)-min(ts(ind1)))<=60 );ind=ind1(ind2);
tsf=seconds(ts(ind)-min(ts(ind)) );tff=seconds(tf(ind)-min(ts(ind)) );n1f=n1(ind);n2f=n2(ind);signsf=signs(ind);durations=(tff-tsf);

%% Converting node labels to numbers
nodes=unique([n1f' n2f']);
n1f_=zeros(length(n1f),1);n2f_=zeros(length(n2f),1);
L1=length(n1f);
countij=0;indij=zeros(10,1);
L2=length(nodes);N=length(nodes);

%% Final variables with no self loops 
[indij, n1f_ , n2f_]=cal_indij(n1f,n2f,L1,L2,nodes,n1f_,n2f_,countij,indij);
start_no_self=tsf(indij);end_no_self=tff(indij);
node1_no_self=n1f_(indij);node2_no_self=n2f_(indij);
signs_no_self=signs(indij);durations_no_self=durations(indij)+1;
signs_no_self=lower(signs_no_self);signs_numeric = iden_sign(signs_no_self);
SHUFFLE=1;
%% MOTIF counting

[motif_counts,  motif_times]...
=count_motifs_deltaT(SHUFFLE,node1_no_self,node2_no_self,start_no_self,end_no_self,durations_no_self,signs_numeric,THRESHOLD_TIME_GAP);
  
if ismember(id,indF)==1
    for mm=1:13
        cc=count_motif_F_null(mm);
        store_motif_times_F_null(mm,cc+1:cc+motif_counts(mm))=motif_times(mm,1:motif_counts(mm));
        count_motif_F_null(mm)=count_motif_F_null(mm)+motif_counts(mm);
    end

else
    for mm=1:13
        cc=count_motif_S_null(mm);
        store_motif_times_S_null(mm,cc+1:cc+motif_counts(mm))=motif_times(mm,1:motif_counts(mm));
        count_motif_S_null(mm)=count_motif_S_null(mm)+motif_counts(mm);
    end

end

SHUFFLE=0;
%% MOTIF counting
[motif_counts,  motif_times]...
=count_motifs_deltaT(SHUFFLE,node1_no_self,node2_no_self,start_no_self,end_no_self,durations_no_self,signs_numeric,THRESHOLD_TIME_GAP);
  
if ismember(id,indF)==1
    for mm=1:13
        cc=count_motif_F(mm);
        store_motif_times_F(mm,cc+1:cc+motif_counts(mm))=motif_times(mm,1:motif_counts(mm));
        count_motif_F(mm)=count_motif_F(mm)+motif_counts(mm);
    end

else
    for mm=1:13
        cc=count_motif_S(mm);
        store_motif_times_S(mm,cc+1:cc+motif_counts(mm))=motif_times(mm,1:motif_counts(mm));
        count_motif_S(mm)=count_motif_S(mm)+motif_counts(mm);
    end

end
end

%%
ind9=[1,3,2,12,13,7,8,9,10];
ind4=[11,4,6,5];

h=figure(4); opq=.4;
x_null=(count_motif_S_null(ind9)+count_motif_F_null(ind9))./...
    (sum((count_motif_S_null(ind9)))+sum((count_motif_F_null(ind9))));
x=(count_motif_S(ind9)+count_motif_S(ind9))/...
    (sum((count_motif_S(ind9)))+sum((count_motif_S(ind9)))); 
bar(x./x_null,'facealpha',opq,'BaseValue',1);

xticks([1:9]); str={'ab-ba','ab-bo','ab-by','ao-xo','ao-xa','ao-xy','ab-xo',...
     'ab-xa','ab-xb','ab-xy','ao-ay','ab-ao','ab-ay'};
  yticks(linspace(0,10,11));
 xticklabels(upper(str(ind9))); 
 ylabel('$\frac{p_{empirical}}{p_{null}}$','interpreter','latex','FontSize',20);
 set(gcf,'Position',[40 40 800 400]);
      title('Pairwise Interactions','interpreter','latex','FontSize',20);box off;
      yline(1,'k');
 print(gcf,'paper_Pshift_pair_REALvsNULL.png','-dpng','-r300');
  
 set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h,'svgs/fig3a_deltaT.svg','-dsvg','-r0');

 %%
 h=figure(5);
 a=5;b=3;
x_null=(count_motif_S_null(ind4)+count_motif_F_null(ind4))./...
    (sum((count_motif_S_null(ind4)))+sum((count_motif_F_null(ind4))));
x=(count_motif_S(ind4)+count_motif_S(ind4))/...
    (sum((count_motif_S(ind4)))+sum((count_motif_S(ind4)))); 
bar(x./x_null,'facealpha',opq,'BaseValue',1);
  xticks([1:4]);  xticklabels(upper(str(ind4))) 
  yticks(linspace(0,10,11));
 ylabel('$\frac{p_{empirical}}{p_{null}}$','interpreter','latex','FontSize',20);
     set(gcf,'Position',[800 40 400 400]);
     title('Group Interactions','interpreter','latex','FontSize',20);box off;
print(gcf,'paper_Pshift_group_REALvsNULL.png','-dpng','-r300');

  set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h,'svgs/fig5a_deltaT.svg','-dsvg','-r0');
% close all
%%
%% functions
function e=TnB(i,j)
e=10*i+j;
end
function e=flipit(e0) %% flips the edge array: source exhancged with target
    e=str2num(fliplr(num2str(e0)));
end



function [indij, n1f_ , n2f_]=cal_indij(n1f,n2f,L1,L2,nodes,n1f_,n2f_,countij,indij) 
%% This function filters the game operators interactions with the players.
for i=1:L1
t1=n1f{i};
t2=n2f{i};

    
    if t1 == "Op" || t2 == "Op" 
    
    else
    for j=1:L2
    t3=nodes{j};
        if strcmp(t1,t3)==1
        n1f_(i)=j;
        end
        if strcmp(t2,t3)==1
        n2f_(i)=j;
        end
    end
        if n2f_(i)~= n1f_(i)
            countij=countij+1;
            indij(countij)=i;
        end
    end
    
end
end
