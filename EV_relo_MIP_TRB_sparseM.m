function EV_relo_MIP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code for MA Tai-Yu, PANTELIDIS Theodoros, CHOW Joseph Y.J. Optimal
% queueing-based rebalancing for one-way electric carsharing systems
% with stochastic demand. Proceeding of the 2019 TRB Annual Meeting.
% Washington D.C: Transportation Research Board (TRB), 2019, 17 p.
% Tai-Yu Ma 28/09/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIP formulation solved by matlab MIP solver
% the problem is reformulated in a graph with N*H nodes, each node is a
% point in (i,g) plan i belongs to N (geo cood), g belongs to H (charging levels) 
% the sigle-commodity flow problem is solved by the link-based formulation
% use sparce matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment setting
% 1 : 6 node case 
% 2 : large network: 
% number of decision variables is NH(NH+m)+|A| 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%parameters
clear all

casestudy=1;
temp=cputime;
%load 'Lamda_400.txt'% 100 nodes and 4 levels 

    
if casestudy==1%case study N=6
   
    %Lamda=[0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 8 2 2    0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 15 0.1 ]';
    
    load Lamda6; %arrival rates at each node-charge i

   N=6%set of nodes
   H=4%charge levels
   O = 3% n of origins 
   F = 3 % num of vehs
   NJ = 2; % num of charging station 
   set_O=[3 7 14]'%origin nodes of vehicles
   set_J=[2 6]'% set of charging station nodes
   yig  =[1 1 1]'
   Cj=3% max num of veh at a node
   rho=[ 0.2236 0.6416 1.1576]'%rho generated based on b=0;alfa=0.95;

   yig_vec = zeros(N*H,1);
   yig_vec(set_O)=1;%n of idle veh at origins
   dist_ij=5%time between two grid neighbor nodes
   dist_h=25% time between two charge neighbor nodes
   activate_queueConstr=0 % 1 means non-myopic model, 0 means myopic model 
   NN1=N*H % size_i is total number of nodes-charge in the node-charge graph
   B=F;
   big_M = 10000;
   Cap = 3;% capacity of charging station : if hetero, change to vec
   theta=0.2;
  
   mu= 1/10*ones(NN1,1);%service rate : 1 customer per 10 minutes

else %large case study 
   N= 100%set of nodes
   H=4%charge levels
   O=10% n of origins 
   F =10 % num of vehs
   %  set_O=[1:10  ]'%for N=10
   %set_O=[1 3 5 7 8 11 12 16 18 20  ]'%for N=20
   set_O=[7 8 16 18 20 27 37 38 45 50 ]'
   NJ = 4; % num of charging station n
   %set_J=[2 5 7 9]'%for N=10
   % set_J=[2 5 17 49]'% set of charging station nodes
   set_J=[15 27 36 49]'% set of charging station nodes
   yig  =ones(F,1);
   yig_vec = zeros(N*H,1);
   yig_vec(set_O)=1;%n of idle veh at origins
   dist_ij=5%time between two grid neighbor nodes
   dist_h=25% time between two charge neighbor nodes
   activate_queueConstr=0
   NN1=N*H % size_i is total nodes-charge in the system
   B=F
   big_M = 10000;
   Cap = 4 ;% capacity of charging station: if hetero, change to vec 
   theta=0.2;
   Cj=3% max num of veh at a node
   rho=[ 0.2236 0.6416 1.1576]'%rho generated based on b=0;alfa=0.95;
   mu= 0.6*N*ones(NN1,1);%service rate : 1 customer per 20 minutes
   load  Lamda100;
end

if O~= size(set_O,1)
    error('O~= size(set_O,1) error ')
end
if NJ~= size(set_J,1)
    error('NJ~= set_J(set_O,1) error ')
end
if H<2
    error('H<2')
end
%%%%%%%%%
%generate graph 
%%%%%%%

n_arc=0;% n of arcs in the node-charge graph 
MAX_INOUTARCS=8; % max num of outgoing arcs of a node in G  
r_ij = zeros(NN1*MAX_INOUTARCS,1);
arcset_in  = zeros(NN1,MAX_INOUTARCS);
arcset_out = zeros(NN1,MAX_INOUTARCS);
arcset_charging = zeros(NN1*2,1);
count_tract_in  = zeros(NN1  ,1);
count_tract_out = zeros(NN1  ,1);
n_arc_charging=0;

II=[];JJ=[];VV=[];

G_node1=[];G_node2=[];%for visualization
for i=1:H
    for j=1:N-1
     index_i = (i-1)*N+j;   
     index_j = (i-1)*N+j+1; 
     G_node1= [G_node1;index_i]; G_node2= [G_node2;index_j];

     %set up link (i,j)
     n_arc=n_arc+1;
     r_ij(n_arc)= dist_ij;
     count_tract_in(index_j)=count_tract_in(index_j)+1;
     an_index1= count_tract_in(index_j);
     arcset_in(index_j,an_index1) = n_arc;% n_arc is arc id
     
     count_tract_out(index_i)=count_tract_out(index_i)+1;
     an_index2= count_tract_out(index_i);
     arcset_out(index_i,an_index2)= n_arc;
     II=[II;index_i ];JJ=[JJ;index_j];VV=[VV;n_arc];
     
     %set up link (j,i)
     G_node1= [G_node1;index_j]; G_node2= [G_node2;index_i];
     n_arc=n_arc+1;
     r_ij(n_arc)= dist_ij;
     count_tract_in(index_i)=count_tract_in(index_i)+1;
     an_index1= count_tract_in(index_i);
     arcset_in(index_i,an_index1) = n_arc;% n_arc is arc id
     
     count_tract_out(index_j)=count_tract_out(index_j)+1;
     an_index2= count_tract_out(index_j);
     arcset_out(index_j,an_index2)= n_arc;
     II=[II;index_j];JJ=[JJ;index_i];VV=[VV;n_arc];
    end
end

for i=1:NJ % set up links on charging stations
    for h=1:H-1
        index_i = set_J(i,1)+(h-1)*N;%start node
        index_j = set_J(i,1)+ h*N;   %end node
        G_node1= [G_node1;index_i]; G_node2= [G_node2;index_j];
        %set up link (i,j)
        n_arc=n_arc+1;
        r_ij(n_arc)= dist_h; %go up
        count_tract_in(index_j)=count_tract_in(index_j)+1;
        an_index1= count_tract_in(index_j);
        arcset_in(index_j,an_index1) = n_arc;% n_arc is arc id
        count_tract_out(index_i)=count_tract_out(index_i)+1;
        an_index2= count_tract_out(index_i);
        arcset_out(index_i,an_index2)= n_arc;
        II=[II;index_i];JJ=[JJ;index_j];VV=[VV;n_arc];
        n_arc_charging = n_arc_charging+1;
        arcset_charging(n_arc_charging)= n_arc;
  
    end
end
n_arc % check
link_dic = sparse(II,JJ,VV,NN1,NN1);% store node-link dictionary

%% checked ok

% generate charging levels 
Level=zeros(NN1,1);
for h=1:H
   h_vec= (h-1)*N+1: (h-1)*N+N;
   Level(h_vec,1)=h;
end

%%%%%%%%%%%%%%%%%%%
%check
%%%%%%%%%%%%%%%%%%%%%%
if size(rho,1)~=Cj
    error(size(rho,1)~=Cj)
end
if F<=0
    error(' F<=0');
end

%compute t_ij % geo distance between nodes 
dist_t=zeros(N,N);
for i=1:N-1
    for j=i+1:N
        dist_t(i,j)=(j-i)*dist_ij;
         dist_t(j,i)=dist_t(i,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dist_tij is dist between nodes-charge
%%%%%%%%%%%%%%%%%%%%%%%%%%
dist_tij= zeros(NN1,NN1);
for i=1:NN1-1
    for j=i+1:NN1
        index_i = rem(i,N);
        index_j = rem(j,N);
        if index_i==0 
            index_i=N;
        end
        if index_j==0 
            index_j=N;
        end
        dist_tij(i,j)= dist_t(index_i,index_j);
        dist_tij(j,i)= dist_tij(i,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up the MIP problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_subpath=  H*(H-1)/2;

length_X= NN1^2 + NN1*Cj+ n_arc+NJ*N_subpath  ;%xij, yjm, wij
init_rowi=zeros(1, length_X);
% eq (1)  objective function
f= init_rowi;

%objecive function
for i=1:NN1
    for j=1:NN1
        index = (i-1)*NN1 +j;
        f(1,index)=dist_tij(i,j)*Lamda(i,1);
    end
end
f(1, NN1*(NN1+Cj)+1 : NN1*(NN1+Cj)+ n_arc) = theta*r_ij(1:n_arc);

%constraints
%Aeq=[];beq=[];A=[];b=[];

N_A = NN1*(Cj-1)+NN1+NN1*NN1+O+2*(NN1-O)+N;%larger than necessary
N_Aeq = 2*NN1+1+NN1+N;

%A=zeros(N_A,length_X);
b=zeros(N_A,1);
%Aeq=zeros(N_Aeg,length_X);
beq=zeros(N_Aeq,1);

Row_Aeq=0;
Row_A=0;
%set up sparse matrix
count_index_A=0;
count_index_Aeq=0;
N_MAX= NN1.*NN1+NN1.*(Cj+5)-O+NN1;
I_A = zeros(N_MAX,1);
J_A = zeros(N_MAX,1);
V_A = zeros(N_MAX,1);
I_Aeq = zeros(N_MAX,1);
J_Aeq = zeros(N_MAX,1);
V_Aeq = zeros(N_MAX,1);

%eq (2)  sum Xij=1 

for i=1:NN1
    %rowi=init_rowi;
    Row_Aeq = Row_Aeq+1;
    for j=1:NN1
        if Level(j)>=Level(i) 
           count_index_Aeq = count_index_Aeq+1;
           index = (i-1)*NN1+j;
         %  rowi(index)=1;
           I_Aeq(count_index_Aeq,1) = Row_Aeq;
           J_Aeq(count_index_Aeq,1) = index;
           V_Aeq(count_index_Aeq,1) = 1;
        end
    end
    %Aeq=[Aeq;rowi];
    %beq=[beq; 1];
    beq(Row_Aeq,1)=1;
end

%eq (3) sum Xij=0

for i=N+1:NN1
  % rowi=init_rowi;
   Row_Aeq = Row_Aeq+1;
   for j=1:NN1
       if Level(j)<Level(i)  
         count_index_Aeq = count_index_Aeq+1;
         index = (i-1)*NN1+j;        
         %   rowi(index)=1;
         I_Aeq(count_index_Aeq,1) = Row_Aeq;
         J_Aeq(count_index_Aeq,1) = index;
         V_Aeq(count_index_Aeq,1) = 1;
        end
    end
    %Aeq=[Aeq;rowi];
  % beq=[beq;0];
    beq(Row_Aeq,1)=0;
end
% test
 
%eq 4  Yjm<Yjm_1 N*m
for j=1:NN1
    for m=2:Cj
        %rowi=init_rowi;
         Row_A = Row_A +1;
         
         count_index_A = count_index_A+1;
         index1= NN1*NN1 + (j-1)*Cj+m-1;
        %rowi(index1)=-1;
         I_A(count_index_A,1) = Row_A;
         J_A(count_index_A,1) = index1;
         V_A(count_index_A,1) = -1;
       
        %rowi(index1+1)=1;
         count_index_A = count_index_A+1;
         I_A(count_index_A,1) = Row_A;
         J_A(count_index_A,1) = index1+1;
         V_A(count_index_A,1) = 1;
        %A=[A;rowi];
        %b=[b; 0];
        b(Row_A,1)=0;
    end
end
%size(A,1) 
%size(b,1)

% eq 5 sum lamda Xij<=mu(Yjm*rho+sum(Yjm*(rho_m-ro_m_1)
if activate_queueConstr==1
   for j=1:NN1
     %rowi=init_rowi;
     Row_A = Row_A +1;
     for i=1:NN1
        index1=(i-1)*NN1+j;
        %rowi(index1)=lamda(i);
        count_index_A = count_index_A+1;
        I_A(count_index_A,1) = Row_A;
        J_A(count_index_A,1) = index1;
        V_A(count_index_A,1) = Lamda(i);
    end
    %m loop
    index2= NN1*NN1+(j-1)*Cj+1;
    %rowi(index2)= -mu(j).*rho(1);
     count_index_A = count_index_A+1;
     I_A(count_index_A,1) = Row_A;
     J_A(count_index_A,1) = index2;
     V_A(count_index_A,1) = -mu(j).*rho(1);
    for m=2:Cj 
        index= NN1*NN1+(j-1)*Cj+m;
       % rowi(index)= mu(j).*(-rho(m)+rho(m-1));
        count_index_A = count_index_A+1;
        I_A(count_index_A,1) = Row_A;
        J_A(count_index_A,1) = index;
        V_A(count_index_A,1) = mu(j).*(-rho(m)+rho(m-1));
    end
    %A=[A;rowi];
    %b=[b;0];
     % A(Row_A,:)=rowi;
      b(Row_A,1)=0;
   end
end    
    
% eq (6) sum Yjm=B

%rowi=init_rowi;
Row_Aeq = Row_Aeq+1;
for j=1:NN1
    for m=1:Cj
        index_x= NN1*NN1+(j-1)*Cj+m;
        %rowi(index_x)=1;
        count_index_Aeq = count_index_Aeq+1;
        I_Aeq(count_index_Aeq,1) = Row_Aeq;
        J_Aeq(count_index_Aeq,1) = index_x;
        V_Aeq(count_index_Aeq,1) = 1;
    end
end
%Aeq=[Aeq;rowi];
%beq=[beq; B];
beq(Row_Aeq,1)=B;


% eq (7)  Xij<=Yjm1
for i=1:NN1
    for j=1:NN1
        %rowi=init_rowi;
         Row_A = Row_A +1;
        index1= (i-1)*NN1+j;
        index2=NN1*NN1+(j-1)*Cj+1;
        %rowi(index1)=1; 
        count_index_A = count_index_A+1;
        I_A(count_index_A,1) = Row_A;
        J_A(count_index_A,1) = index1;
        V_A(count_index_A,1) = 1;
        %rowi(index2)=-1;
        count_index_A = count_index_A+1;
        I_A(count_index_A,1) = Row_A;
        J_A(count_index_A,1) = index2;
        V_A(count_index_A,1) = -1;
        %A=[A;rowi];
        %b=[b; 0];
         b(Row_A,1)=0;
    end
end
%size(A,1) 
%size(b,1)

% eq (8) sum -inflow + sum outflow <= yjm*bigM
ss=1:NN1;
ss(set_O)=[];
for i=1:size(ss,2)
    %rowi=init_rowi;
    Row_A = Row_A +1;
    nodeid= ss(i);
    %incoming flow
    for j=1:MAX_INOUTARCS
        if arcset_in(nodeid,j) > 0
           linkid= arcset_in(nodeid,j);
           index = NN1*NN1+NN1*Cj+linkid;
           %rowi(index)= 1;
           count_index_A = count_index_A+1;
           I_A(count_index_A,1) = Row_A;
           J_A(count_index_A,1) = index;
           V_A(count_index_A,1) = 1;
        else
            break;
        end
    end
     %outgoing flow
     for j=1:MAX_INOUTARCS
        if arcset_out(nodeid,j) > 0
           linkid= arcset_out(nodeid,j);
           index = NN1*NN1+NN1*Cj+linkid;
           %rowi(index)= -1;
           count_index_A = count_index_A+1;
           I_A(count_index_A,1) = Row_A;
           J_A(count_index_A,1) = index;
           V_A(count_index_A,1) = -1;
        else
            break;
        end
     end    
    index2= NN1*NN1+ (nodeid-1)*Cj+1;
    %rowi(index2) = - big_M ;%big_M
    count_index_A = count_index_A+1;
    I_A(count_index_A,1) =  Row_A  ;
    J_A(count_index_A,1) =  index2 ;
    V_A(count_index_A,1) = -big_M  ;
   
    %A=[A;rowi];
    %b=[b;  0];
    b(Row_A,1)=0;
end


%eq (9) - (sum inflow - sum outflow) < = yjm*bigM
for i=1:size(ss,2)
   % rowi=init_rowi;
    Row_A = Row_A +1;
    nodeid= ss(i);
    %incoming flow
    for j=1:MAX_INOUTARCS
        if arcset_in(nodeid,j) > 0
           linkid= arcset_in(nodeid,j);
           index = NN1*NN1+NN1*Cj+linkid;
           %rowi(index)= -1;
           count_index_A = count_index_A+1;
           I_A(count_index_A,1) = Row_A;
           J_A(count_index_A,1) = index;
           V_A(count_index_A,1) = -1;
        else
            break;
        end
    end
     %outgoing flow
     for j=1:MAX_INOUTARCS
        if arcset_out(nodeid,j) > 0
           linkid= arcset_out(nodeid,j);
           index = NN1*NN1+NN1*Cj+linkid;
          % rowi(index)= 1;
           count_index_A = count_index_A+1;
           I_A(count_index_A,1) = Row_A;
           J_A(count_index_A,1) = index;
           V_A(count_index_A,1) = 1;
        else
            break;
        end
     end    
    index2= NN1*NN1+ (nodeid-1)*Cj+1;
    %rowi(index2) = - big_M;%big_M
    count_index_A = count_index_A+1;
    I_A(count_index_A,1) =  Row_A  ;
    J_A(count_index_A,1) =  index2 ;
    V_A(count_index_A,1) = -big_M  ;
    %A=[A;rowi];
    %b=[b;  0];
    b(Row_A,1)=0;
end
%size(A,1) 
%size(b,1)

% eq (10) sum Wij = sum Yjm
for j=1:NN1
    %rowi=init_rowi;
    Row_Aeq = Row_Aeq+1;
    %incoming flow
    for i=1:MAX_INOUTARCS
       if arcset_in(j,i) > 0
           linkid= arcset_in(j,i);
           index = NN1*NN1+NN1*Cj+linkid;
          % rowi(index)= 1;
           count_index_Aeq = count_index_Aeq+1;
           I_Aeq(count_index_Aeq,1) = Row_Aeq;
           J_Aeq(count_index_Aeq,1) = index;
           V_Aeq(count_index_Aeq,1) = 1;
        else
            break;
        end
    end
     %outgoing flow
    for i=1:MAX_INOUTARCS
        if arcset_out(j,i) > 0
           linkid= arcset_out(j,i);
           index = NN1*NN1+NN1*Cj+linkid;
          %rowi(index)= -1;
           count_index_Aeq = count_index_Aeq+1;
           I_Aeq(count_index_Aeq,1) = Row_Aeq;
           J_Aeq(count_index_Aeq,1) = index;
           V_Aeq(count_index_Aeq,1) = -1;
        else
            break;
        end
    end 
        
    for m=1:Cj
        index = NN1*NN1+ (j-1)*Cj+m;
        %rowi(index)= -1;
         count_index_Aeq = count_index_Aeq+1;
         I_Aeq(count_index_Aeq,1) = Row_Aeq;
         J_Aeq(count_index_Aeq,1) = index;
         V_Aeq(count_index_Aeq,1) = -1;
    end
    %Aeq=[Aeq;rowi];
    %beq=[beq; -yig_vec(j) ];
    beq(Row_Aeq,1)= -yig_vec(j);
end
%size(Aeq,1) 
%size(beq,1)

% eq (11) sum prs =wij
id_arc_charge=0;
for i=1:NJ % set up links on charging stations
    for arc_ch=1:H-1 % arc
        Row_Aeq = Row_Aeq+1;
        index_p=0;
        g=arc_ch;
        h=arc_ch+1;
        id_arc_charge=id_arc_charge+1;
        for gg=1:H-1
            for hh=gg+1:H
                index_p = index_p+1;
                if gg<=g && hh>=h
                   index_i = NN1^2 + NN1*Cj + n_arc + (i-1)*N_subpath+index_p; %N_subpath = H(H-1)/2
                   count_index_Aeq = count_index_Aeq+1;
                   I_Aeq(count_index_Aeq,1) = Row_Aeq;
                   J_Aeq(count_index_Aeq,1) = index_i;
                   V_Aeq(count_index_Aeq,1) = 1;
                end
            end
        end
        link_id = arcset_charging(id_arc_charge);
        index_w = NN1^2 + NN1*Cj+ link_id;
        count_index_Aeq = count_index_Aeq+1;
        I_Aeq(count_index_Aeq,1) = Row_Aeq;
        J_Aeq(count_index_Aeq,1) = index_w;
        V_Aeq(count_index_Aeq,1) = -1;
        beq(Row_Aeq,1)=0;
    end
end

% eq (12) % sum pij <=uj 
for i=1:NJ % 
    Row_A = Row_A+1;
    index_p=0;        
    for g=1:H-1
        for h=g+1:H
            index_p =index_p+1;
            index_i = NN1^2 + NN1*Cj + n_arc + (i-1)*N_subpath+index_p; 
            count_index_A = count_index_A+1;
            I_A(count_index_A,1) = Row_A;
            J_A(count_index_A,1) = index_i;
            V_A(count_index_A,1) = 1;
        end
    end
    b(Row_A,1)=Cap;
end


intcon =[1:length_X ];
lb = zeros(length_X,1);
ub = [ones(NN1*(NN1+Cj),1); Inf*ones(length_X-NN1*(NN1+Cj),1)];

I_A=I_A(1:count_index_A,1);
J_A=J_A(1:count_index_A,1);
V_A=V_A(1:count_index_A,1);

I_Aeq=I_Aeq(1:count_index_Aeq,1);
J_Aeq=J_Aeq(1:count_index_Aeq,1);
V_Aeq=V_Aeq(1:count_index_Aeq,1);

A   = sparse(I_A,J_A,V_A,      Row_A, length_X);
Aeq = sparse(I_Aeq,J_Aeq,V_Aeq, Row_Aeq, length_X);
b=b(1:Row_A,1);
beq=beq(1:Row_Aeq , 1);

t=cputime;
disp('eq matrix preparation')
diff =t -temp
temp=t;
length_X %check

size(A,1)+size(Aeq,1)%check

options = optimoptions(@intlinprog,'MaxNodes',10000)
[x,fval,exitflag,output] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options);

t=cputime;
disp('exec time')
diff =t -temp
fval

% get solution
Xij=zeros(NN1,NN1);%
Yjm=zeros(NN1,Cj);
Wij=zeros(n_arc,1);%node based, sparse matrix
Pij=zeros(NJ*N_subpath,1);% subpath flow
for i=1:NN1
    for j=1:NN1
        index= (i-1)*NN1+j;
        Xij(i,j)=x(index);
    end
end
Xij;

for j=1:NN1
    for m=1:Cj
        index= NN1*NN1+(j-1)*Cj+m;
        Yjm(j,m) = x(index);
    end
end
Yjm ;

for i=1:n_arc
    index = NN1*NN1+NN1*Cj+i;
    Wij(i)=x(index);
end


for i=1:NJ*N_subpath
    index = NN1*NN1+NN1*Cj+n_arc+i;
    Pij(i,1)=x(index);
end

Wij;
z=f*x
Yjm ;
Pij


% visu flow 
if  casestudy==1  % 6 node case
    figure (1)

    names = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24'};
    G = digraph(G_node1,G_node2,Wij',names)
    LWidths = 5*G.Edges.Weight+1/max(G.Edges.Weight);
    p = plot(G,'Layout','force','EdgeLabel',G.Edges.Weight, 'LineWidth',LWidths)

    xdata=[];ydata=[];
    for h=1:4
        for j=1:N
            cood_x =  (j);xdata=[xdata;cood_x];
            cood_y = (h);ydata=[ydata;cood_y];
        end
    end

    p.XData= xdata';
    p.YData= ydata';
    p.NodeColor = 'r'
    p.MarkerSize = 8;
    grid on
    xlabel('Node (i)')
    ylabel('Charge level (h)')

end
  

%  clear s
 % s = int2str(1);
 %   for i=2:NN1
  %      chr = int2str(i)
  %      aa=','
  %      s = strcat(s,aa,chr)
 %   end

%check in-out arc sets

%TEST in-/out- arc sets 
%for i=1:24
  %  aa1=[];
 %   for j=1:5
 %       if arcset_out(i,j)>0
  %          aa1 = [aa1;arcset_out(i,j)];
  %      else
   %         break;
   %     end
  %  end
 %   subplot(4,6,i);
 %   G = digraph(G_node1(aa1),G_node2(aa1));
% %   p = plot(G)   
  %  clear G
%end
         




