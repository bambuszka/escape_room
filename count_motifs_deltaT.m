
function [motif_counts , motif_times]=...
count_motifs_deltaT(SHUFFLE,node1_no_self,node2_no_self,start_no_self,end_no_self,durations_no_self,signs_numeric,THRESHOLD_TIME_GAP)
motif_times=zeros(13,1);motif_counts=zeros(13,1);
[elist_new,Tstart_new,Tend_new]=...
create_elist_new(node1_no_self,node2_no_self,start_no_self,end_no_self,durations_no_self,signs_numeric);
n_events=length(Tstart_new)
% delta_T_between_succesive_interactions=Tstart_new(2:n_events,:)-Tend_new(1:n_events-1,:);
delta_T_between_succesive_interactions=Tstart_new(2:n_events,:)-Tstart_new(1:n_events-1,:);

Enew=length(elist_new);
%1 abba,  2  abbo, 3 abby, 4 aoxo, 5 aoxa,%6 aoxy,
%7 abxo,%8 abxa,%9 abxb,%10 abxy,%11 aoay,%12 abao,%13 abay
if SHUFFLE==1;shuffled=randperm(Enew);
else;shuffled=[1:Enew];end
for i=1:Enew-1
    e=(shuffled(i)); if (e<Enew) & (delta_T_between_succesive_interactions(i)<=THRESHOLD_TIME_GAP);
    a=elist_new(e,1);b=elist_new(e,2);
    c=elist_new(shuffled(i+1),1);d=elist_new(shuffled(i+1),2);
    
    if b==10
          if c~=a && c~=10 
                if d==10
                    m=4;motif_counts(m)=motif_counts(m)+1;motif_times(m,motif_counts(m))= Tstart_new(e);
                elseif d==a
                    m=5;motif_counts(m)=motif_counts(m)+1;motif_times(m,motif_counts(m))= Tstart_new(e);
                else
                    m=6;motif_counts(m)=motif_counts(m)+1;motif_times(m,motif_counts(m))= Tstart_new(e);
                end
          elseif c==a && d~=a && d~=10
                    m=11;motif_counts(m)=motif_counts(m)+1;motif_times(m,motif_counts(m))= Tstart_new(e);

          end
    
       
    else
    
        if c==b 
                if d==a
                    m=1;motif_counts(m)=motif_counts(m)+1;motif_times(m,motif_counts(m))= Tstart_new(e);
                elseif d==10
                    m=2;motif_counts(m)=motif_counts(m)+1;motif_times(m,motif_counts(m))= Tstart_new(e);
                else
                    m=3;motif_counts(m)=motif_counts(m)+1;motif_times(m,motif_counts(m))= Tstart_new(e);
                end

        elseif c==a
                 if d~=b && d~=10
                    m=13;motif_counts(m)=motif_counts(m)+1;motif_times(m,motif_counts(m))= Tstart_new(e);
                 elseif d==10
                    m=12;motif_counts(m)=motif_counts(m)+1;motif_times(m,motif_counts(m))= Tstart_new(e);
                 end
  %1 abba,  2  abbo, 3 abby, 4 aoxo, 5 aoxa,%6 aoxy,
%7 abxo,%8 abxa,%9 abxb,%10 abxy,%11 aoay,%12 abao,%13 abay               
        elseif c~=a && c~=b
                   if d==a
                        m=8;motif_counts(m)=motif_counts(m)+1;motif_times(m,motif_counts(m))= Tstart_new(e);
                   elseif d==b
                        m=9;motif_counts(m)=motif_counts(m)+1;motif_times(m,motif_counts(m))= Tstart_new(e);
                   elseif d~=10
                       m=10;motif_counts(m)=motif_counts(m)+1;motif_times(m,motif_counts(m))= Tstart_new(e);
                   elseif d==10
                       m=7;motif_counts(m)=motif_counts(m)+1;motif_times(m,motif_counts(m))= Tstart_new(e);
                   end 
        end
    end
end
%1 abba,  2  abbo, 3 abby, 4 aoxo, 5 aoxa,%6 aoxy,
%7 abxo,%8 abxa,%9 abxb,%10 abxy,%11 aoay,%12 abao,%13 abay
end
end