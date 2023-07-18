%% Description and download data

% The Ecoli_case_study.m file demonstrates the case study of the Ecoli
% biomass growth preference.
% data reference: Orth, J. D.; Fleming, R. M. T.; Palsson, B. Ø. Reconstruction and Use of Microbial Metabolic Networks: The Core Escherichia Coli Metabolic Model as an Educational Guide. EcoSal Plus 2010, 4 (1), ecosalplus.10.2.1. https://doi.org/10.1128/ecosalplus.10.2.1.
% Download ECOLI_example.mat file and load data
load ECOLI_example
modelJoint=model;
S=modelJoint.S;
S=full(S);
[m n]=size(S);

% Produce reaction adjacency matrix
S_hat=double(boolean(S ~= 0));
A=transpose(S_hat)*(S_hat);

% check  reactions' irreversible
[m n]=size(S);
irre_check=zeros(n,1);
for i=1:n
    if modelJoint.lb(i,1) < -0.001 && modelJoint.ub(i,1)  > 0.001
        irre_check(i)=1;     
    else
        irre_check(i,1)=0; 
    end
end


%  S2m, the unfolded version of the stoichiometric matrix of the 2m, the forward and reverse reactions.
m=size(A,2);
S_2m=[S, -S]*[eye(m), zeros(m);  zeros(m),diag(irre_check)];
S_2m_pos=0.5*(abs(S_2m)+S_2m);
S_2m_neg=0.5*(abs(S_2m)-S_2m);


% mass flux
[FBAV,~,~] = xlsread('DATA SET S1.xlsx','Core E. coli Model','C2:C191');
% flux in Gluc: [FBAV,~,~] = xlsread('DATA SET S1.xlsx','Core E. coli Model','C2:C191');
% flux in Gluc_lim: [FBAV,~,~] = xlsread('DATA SET S1.xlsx','Core E. coli Model','D2:D191');
v=FBAV;
v_2m = v; % depends on the FBA data

% calculate the FBA solution based weight
j_v=S_2m_pos*v_2m;
Jv=diag(j_v);
V=diag(v_2m);
MFG=transpose(S_2m_pos*V)* pinv(Jv)*(S_2m_neg*V);


% filter the weights in the flux-dependent graph
[m n]=find(MFG ~= 0);
WW=[];
for i=1:size(m,1)
    weight=MFG(m(i),n(i));
    WW=[WW;weight];
end

s=m;
t=n;
for i=1:size(s,1)
    if s(i,:)>95
        s(i,:)=s(i,:)-95;
    else
        s(i,:)=s(i,:);
    end
end

for i=1:size(t,1)
    if t(i,:)>95
        t(i,:)=t(i,:)-95;
    else
        t(i,:)=t(i,:);
    end
end

% convert the weight data to the corresponding rxns direction 
for i=1:size(t,1)
    if t(i,:)>95 && s(i,:)>95
        WW(i,:)=-WW(i,:);
    elseif t(i,:)>95 && s(i,:)<=95
        WW(i,:)=-WW(i,:);
    elseif t(i,:)<=95 && s(i,:)>95
        WW(i,:)=-WW(i,:);
    else
        WW(i,:)=WW(i,:);
    end
end


% plot the flux-dependent graph
S=string(s)';
T=string(t)';
weights = WW';
G = digraph(S,T,weights');
H = plot(G);

%% applying Min-Cut algorithm on the flux-dependent graph

% determine the source and sink nodes and find the Min-Cut
% s:glucose => t:biomass (Gluc case)
weights = WW';
G = digraph(S,T,weights');
H = plot(G);
[mf,GF,cs,ct] = maxflow(G,'28','13');
H.EdgeLabel = {};
highlight(H,GF,'EdgeColor','b','LineWidth',2);

% label the Min-Cut/ or Max flow pathway
highlight(H,'28','NodeColor','b','MarkerSize',10)
highlight(H,'23','NodeColor','g','MarkerSize',10)

% calculate the total sum of the weight for the extracted pathway
sum(table2array(GF.Edges(:,2))) 

% sum of the weight forGLUC case: 53.9025
% sum of the weight forLIMITED GLUC case:  81.6555




