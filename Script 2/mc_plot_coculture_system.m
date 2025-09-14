% The m file demonstrates how to construct a flux-dependent graph with
% specific sources and sinks.


B=[];
I=[];
A=[];
E=[];
% choosing time point to generate m-c graph, as the representation of the weights within the time interval
t=[20 100 119 150 193 250 303 313 323];

for tt=t
    % creating empty matrixes
glucose_edgepath_sum=zeros(4,3);
glucose_edgepath_num=zeros(4,3);
fru_edgepath_sum=zeros(4,4);
fru_edgepath_num=zeros(4,4);

% read the stoichiometric matrix, Sv=0
for i=1:4
load model.mat
S=full(modelJoint.S);

load FBA_V.mat
% FBA_V: the FBA solutions for the case. row: v; column: solution chaning
% though the time
% change the column # of FBA_V, here, we use time point = [20 100 119 150 193 250 303 313 323] as an
% example
FBAV_time=FBA_V(:,tt);  

% fix biomass, ATP, glucose and fructose uptake
fixindex=[1100 1702 3369 3978 4668 1476 4884 3752 4639 4699 4855 4917 4557 4695 4773 4913 5046 5041]';

%update the bounds for  biomass, ATP, glucose and fructose uptake
updatebound=FBAV_time(fixindex);
modelJoint.lb(fixindex)=updatebound;


rxnsList = { 'ExCom_btoh[u]' 'ExCom_2ppoh[u]' 'ExCom_ac[u]'  'ExCom_etoh[u]'};
%    ExCom_btoh[u]   'ExCom_2ppoh[u]'   'ExCom_ac[u]'
% 'ExCom_acetone[u]'  {'ExCom_etoh[u]'}
obj_v=[];


% build the flux-dependent graph
    modelJoint = changeObjective (modelJoint, rxnsList{i}); 
    FBA = optimizeCbModel (modelJoint, 'max');
    V_obj=FBA.v;
    obj_v=[obj_v, V_obj];


    % check irreversible
    irre_check=modelJoint.rev;

    S=full(modelJoint.S);
    m=size(S,2);
    S_2m=[S, -1.*S]*[eye(m), zeros(m);  zeros(m),diag(irre_check)];
    S_2m_pos=0.5*(abs(S_2m)+S_2m);
    S_2m_neg=0.5*(abs(S_2m)-S_2m);

    w_pos=sum(S_2m_pos,2);
    w_neg=sum(S_2m_neg,2);

    v=obj_v(:,1);
    v_2m = 0.5.*[abs(v)+v ; abs(v)-v];
    j_v=S_2m_pos*v_2m;
    Jv=diag(j_v);
    V=diag(v_2m);
    MFG=transpose(S_2m_pos*V)* pinv(Jv)*(S_2m_neg*V);

% label the nodes (reactions) with number
% filter MFG
    [m n]=find(MFG ~= 0);
% [m n]=find(MFG >= 1);
% filter weight
    WW=[];
    for j=1:size(m,1)
         weight=MFG(m(j),n(j));
         WW=[WW;weight];
    end
    s=m;
    t=n;
    for j=1:size(s,1)
        if s(j,:)>5141
             s(j,:)=s(j,:)-5141;
        else
             s(j,:)=s(j,:);
        end
    end

    for j=1:size(t,1)
        if t(j,:)>5141
             t(j,:)=t(j,:)-5141;
         else
             t(j,:)=t(j,:);
         end
    end

    S=string(s)';
    T=string(t)';
    weights = WW';
    G = digraph(S,T,weights');


% the number denotes as the index of the reaction. e.g., 4985 --> ipa;  5012--> butyrate 
    rxnsnode = { '5011' '4985' '4987'  '5034'};
% rxnsList = { 'ExCom_btoh[u]' 'ExCom_2ppoh[u]' 'ExCom_ac[u]' 'ExCom_etoh[u]'};

figure(i);
H = plot(G);

% source: glucose
% sink: 4639 --> cac pure
[mf,GF,cs,ct] = maxflow(G,'4639',rxnsnode{i} ); 
H.EdgeLabel = {};
highlight(H,GF,'EdgeColor','b','LineWidth',2);
% st = GF.Edges.EndNodes;
% labeledge(H,st(:,1),st(:,2),GF.Edges.Weight);
highlight(H,'4639','NodeColor','b','MarkerSize',10)
highlight(H,rxnsnode{i},'NodeColor','g','MarkerSize',10)
cac_glucose=sum(table2array(GF.Edges(:,2))) ;
cac_glucose_num=size(GF.Edges,1);
% sink: 3362 --> cacfused 
[mf,GF,cs,ct] = maxflow(G,'4855',rxnsnode{i}); 
highlight(H,GF,'EdgeColor','r','LineWidth',2);
highlight(H,'4855','NodeColor','r','MarkerSize',10)
cacfused_glucose= sum(table2array(GF.Edges(:,2)));
cacfused_glucose_num=size(GF.Edges,1);
% sink: 4159 ---> clj fused
[mf,GF,cs,ct] = maxflow(G,'4917',rxnsnode{i});
highlight(H,GF,'EdgeColor','k','LineWidth',2);
highlight(H,'4917','NodeColor','k','MarkerSize',10)
cljfused_glucose= sum(table2array(GF.Edges(:,2))) ;
cljfused_glucose_num=size(GF.Edges,1);

% collect the edge number and sum of the weights
glucose_edgepath_sum(i,:) =[cac_glucose,cacfused_glucose, cljfused_glucose];
glucose_edgepath_num(i,:) =[cac_glucose_num,cacfused_glucose_num, cljfused_glucose_num];


% source: fructose
% sink: 4557 --> cac pure
figure(4+i);
H = plot(G);
[mf,GF,cs,ct] = maxflow(G,'4557',rxnsnode{i} ); 
H.EdgeLabel = {};
highlight(H,GF,'EdgeColor','b','LineWidth',2);
% st = GF.Edges.EndNodes;
% labeledge(H,st(:,1),st(:,2),GF.Edges.Weight);
highlight(H,'4557','NodeColor','b','MarkerSize',10)
highlight(H,rxnsnode{i},'NodeColor','g','MarkerSize',10)
cac_fru=sum(table2array(GF.Edges(:,2))) ;
cac_fru_num=size(GF.Edges,1);
% sink: 4773 --> cacfused 
[mf,GF,cs,ct] = maxflow(G,'4773',rxnsnode{i});
highlight(H,GF,'EdgeColor','r','LineWidth',2);
highlight(H,'4773','NodeColor','r','MarkerSize',10)
cacfused_fru= sum(table2array(GF.Edges(:,2)));
cacfused_fru_num=size(GF.Edges,1);
% sink: 4913 ---> clj fused
[mf,GF,cs,ct] = maxflow(G,'4913',rxnsnode{i});
highlight(H,GF,'EdgeColor','k','LineWidth',2);
highlight(H,'4913','NodeColor','k','MarkerSize',10)
cljfused_fru= sum(table2array(GF.Edges(:,2))) ;
cljfused_fru_num=size(GF.Edges,1);
% sink: 4695 ---> clj pure
[mf,GF,cs,ct] = maxflow(G,'4695',rxnsnode{i});
highlight(H,GF,'EdgeColor','y','LineWidth',2);
highlight(H,'4695','NodeColor','y','MarkerSize',10)
cljpure_fru= sum(table2array(GF.Edges(:,2))) ;
cljpure_fru_num=size(GF.Edges,1);

% collect the edge number and sum of the weights
fru_edgepath_sum(i,:) =[cac_fru,cacfused_fru, cljfused_fru, cljpure_fru];
fru_edgepath_num(i,:) =[cac_fru_num,cacfused_fru_num, cljfused_fru_num, cljpure_fru_num];
gluc_ed = glucose_edgepath_sum./glucose_edgepath_num;
fru_ed = fru_edgepath_sum./fru_edgepath_num;


buoh = sum(gluc_ed(1,:))+sum(fru_ed(1,:));
ipa = sum(gluc_ed(2,:))+sum(fru_ed(2,:));
ac = sum(gluc_ed(3,:))+sum(fru_ed(3,:));
etoh = sum(gluc_ed(4,:))+sum(fru_ed(4,:));



end

% edge density 
B=[B,buoh]; % butanol
I=[I,ipa]; % IPA
A=[A,ac]; %Acetate
E=[E,etoh]; %Ethanol

end

