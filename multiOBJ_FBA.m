%% Description, data and tool download section

% The multi objective FBA was revised from the DMMM model 
% reference:Foster, C., Charubin, K., Papoutsakis, E. T., & Maranas, C. D. (2021). Modeling Growth kinetics, interspecies cell fusion, and metabolism of a Clostridium acetobutylicum/Clostridium ljungdahlii syntrophic coculture. Msystems, 6(1), e01325-20.
% Download the Cac_Clj_MuliFBA_example.m file
% Download the cobratoolbox
% reference: https://github.com/opencobra/cobratoolbox
% Add path to cobratoolbox
% initCobraToolbox

% load experinmental data
load Cac_Clj_MuliFBA_example

% Cac_Clj_MuliFBA_example.m contains the information of: 
% CACOD and CACtime: Cac OD600 Data
% CLJOD and Clj CLJtime: Clj OD600 Data
% COCOD and COCtime: Co-culture OD600 Data
% CLJ_plot: non-hybrid Clj relative abundnace from flow cytometry experiments
% CAC_plot:  non-hybrid Cac relative abundance from flow cytometry experiments
% MIXED_plot and t_RNA: hybrid cell relative abundance from flow cytometry experiments
% GC_ratio: Cac/Clj ratio calculated from genome copy number
% GC_stDEV and GC_time: Cac/Clj ratio standard deviation
% GC_Cac and GC_Clj genome copy number for fusion parameter conersion to cell/L
% concentration: concentration data

%% data pretreatment section

%Calculate relative abundances and standard deviations  from genome copy 
%number data
Cac_GC = (GC_ratio )./(GC_ratio + 1);
GC_stDEV = GC_stDEV./(GC_ratio + 1);

%interpolate OD600 values
tq = [0:0.1:100];
yq1 = interp1(CACtime,CACOD,tq,'pchip');
yq2 = interp1(CLJtime,CLJOD,tq,'pchip');
yq3 = interp1(COCtime,COCOD,tq,'pchip');


%Calculate time-dependent C. acetobutylicum specific growth rate
for i = 2:length(yq1)
    mu1(i,1) = log(yq1(i)/yq1(i-1))/(tq(i)-tq(i-1));
end
mu1 = min(mu1,3);
k = 2;
for i = 1:201
    muset1_init(i) = min(mu1(k),3.0);
    k = k + 1;
end
%Calculate time-dependent C. ljungdahlii specific growth rate
for i = 2:length(yq2)
    mu2(i,1) = log(yq2(i)/yq2(i-1))/(tq(i)-tq(i-1));
end
k = 2;
for i = 1:201
    muset2_init(i) = min(mu2(k),3.0);
    k = k + 1;
end
%Calculate time-dependent co-culture specific growth rate
for i = 2:length(yq2)
    mu3(i,1) = log(yq3(i)/yq3(i-1))/(tq(i)-tq(i-1));
end
k = 2;
for i = 1:201
    muset3_init(i) = min(mu3(k),3.0);
    k = k + 1;
end

%Initialize storage variable
fullx = zeros(0,12);

%Define fusion parameter
f = 0.00564; 
%Convert fusion parameter to cell/L, assume genome copy number
%concentration is equivalent to cell concentration
muL_conversion = 1e-6;
f_converted = f*100/((GC_Cac/muL_conversion)+(GC_Clj/muL_conversion));

% empty matrix for storing the lower and upper bound
SLB=[]; SUB=[];
% empty matrix for storing the growth rates
MU=[]; MUU=[];
for i = 1:631
    if i == 1
        %initial cell abundance conditions, growth rates
        x0 = [10; 90; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
        tspan = [0 i*0.1];
        mu = [muset1_init(i) muset2_init(i) muset3_init(i)];
        MUU=[MUU;mu];
    else
        %iteratively updated cell abundances, growth rates, integration 
        %over 0.1 hr intervals
        tspan = [(i-1)*0.1 i*0.1];
        x0 = xa(end,:);
        mu = [mu1(i) mu2(i) mu3(i)];
        MUU=[MUU;mu];
    end
    
    %Cell death is not included in growth model
    mu = max(mu,0);
    MU=[MU;mu];
    
    %System of ODEs representing hybrid/non-hybrid cell growth and cell
    %fusion
    g=@(t,x) [mu(1)*max(x(1),0) - f*x(1)*max(x(2),0) + mu(1)*max(x(11),0);
        mu(2)*max(x(2),0) - f*x(1)*max(x(2),0) + mu(2)*x(12);...
    f*x(1)*max(x(2),0) - max(mu(1),0)*max(x(3),0);
    f*x(1)*max(x(2),0) - mu(2)*x(4);
    mu(1)*(x(3)-max(x(5),0));
    mu(2)*(x(4)-x(6));...
    mu(1)*(max(x(5),0)-x(7));
    mu(2)*(x(6)-x(8));
    mu(1)*(x(7)-max(x(9),0));
    mu(2)*(x(8)-x(10));
    mu(1)*(max(x(9),0)-x(11));
    mu(2)*(x(10)-x(12))];
    
    %ODE integration
    [t,xa] = ode45(g,tspan,x0);
    
    %Store solution
    fullx = [fullx;xa(1,:)];
    
    %store ODE timesteps
    if i == 1
        fullt = t(end,:);
    else
        fullt = [fullt;t(end,:)];
    end
end
fullt=[0;fullt];
fullt(632,:)=[];
% for i = 1:50
%     store_ind = find(fullt == i);
%     store_time(i,1) = store_ind(1);
% end
% store_time = flipud(store_time);
% for i = 1:length(store_time)
%     fullt(store_time(i)) = [];
%     fullx(store_time(i),:) = [];
% end

pullt = unique(fullt);
for i = 1:length(pullt)
    tpull = find(pullt(i) == fullt);
    pullx(i,:) = fullx(tpull(1),:);
end
for i = 1:size(fullx,2)
    XX(:,i) = interp1(pullt,pullx(:,i),tq,'pchip');
end
for i = 1:size(XX,1)
    X(i,:) = XX(i,:)./sum(XX(i,:));
end
X(:,13) = X(:,2);
X(:,2) = [];
XX(:,13) = XX(:,2);
XX(:,2) = [];


%Interpolate concentration data
tq = [0:0.1:100];
yq1 = interp1(concentration(:,1),concentration(:,2),tq,'pchip');
yq2 = interp1(concentration(:,1),concentration(:,3),tq,'pchip');
yq3 = interp1(concentration(:,1),concentration(:,4),tq,'pchip');
yq4 = interp1(concentration(:,1),concentration(:,5),tq,'pchip');
yq5 = interp1(concentration(:,1),concentration(:,6),tq,'pchip');
yq6 = interp1(concentration(:,1),concentration(:,7),tq,'pchip');
yq7 = interp1(concentration(:,1),concentration(:,8),tq,'pchip');
yq8 = interp1(concentration(:,1),concentration(:,9),tq,'pchip');
yq9 = interp1(concentration(:,1),concentration(:,10),tq,'pchip');


pp_pchip_gluc = pchip(tq,yq1);
pp_pchip_fruc = pchip(tq,yq2);
%pp_pchip_lac = pchip(tq,yq3);
pp_pchip_ac = pchip(tq,yq3);
%pp_pchip_bdo = pchip(tq,yq5);
pp_pchip_acetoin = pchip(tq,yq4);
pp_pchip_etoh = pchip(tq,yq5);
pp_pchip_isoprop = pchip(tq,yq6);  %%%%%%
pp_pchip_but = pchip(tq,yq7);
pp_pchip_acetone = pchip(tq,yq8);
pp_pchip_btoh = pchip(tq,yq9);

%load acetatebound.mat
data(:,1) = tq(1:1000);
data(:,2) = pp_pchip_gluc.coefs(:,3);
data(:,3) = pp_pchip_fruc.coefs(:,3);
% data(:,4) = pp_pchip_lac.coefs(:,3);
data(:,4) = pp_pchip_ac.coefs(:,3);  
data(:,5) =  pp_pchip_isoprop.coefs(:,3);
data(:,6) = pp_pchip_acetoin.coefs(:,3);
% data(:,8) = pp_pchip_bdo.coefs(:,3);
data(:,7) = pp_pchip_btoh.coefs(:,3);
data(:,8) = pp_pchip_but.coefs(:,3);
data(:,9) = pp_pchip_etoh.coefs(:,3);
data(:,10) = pp_pchip_acetone.coefs(:,3);

%Indecies used to partition glucose (1) and fructose (2) uptake according
%to each stoichiometric model in DMMM framework
Cacpureind = [1 2];
Cljpureind = [2];
Cacfusedind = [1 2];
Cljfusedind = [1 2];

%List of cell types in simulation
species = {'Cacpure';'Cacfused';'Cljfused';'Cacfused';'Cljfused';'Cacfused';'Cljfused';'Cacfused';'Cljfused';'Cacfused';'Cljfused';'Cljpure'};

%0.1 hour timesteps
dt = 0.1;

%Define initial concentrations of each extracellular carbon source/sink
C_gluc0 = concentration(1,2);
C_fruc0 = concentration(1,3);
% C_lac0 = concentration(1,4);
C_ac0 = concentration(1,4);
% C_bdo0 = concentration(1,6);
C_acetoin0 = concentration(1,5);
C_etoh0 = concentration(1,6);
C_isoprop0 = concentration(1,7);
C_but0 = concentration(1,8);
C_acetone0 = concentration(1,9);
C_btoh0 = concentration(1,10);

%Define reaction list, index, and molecular weight for FVA, maximum titer
%calculations
rxnsList = {'ExCom_ac[u]';'ExCom_2ppoh[u]';'ExCom_actn__R[u]';'ExCom_btoh[u]';'ExCom_but[u]';'ExCom_etoh[u]';'ExCom_acetone[u]'};
rxn_index = [4987; 4985; 4991; 5011; 5012; 5034; 4988];
MW = [0.060052; 0.0601; 0.08811; 0.074121; 0.08811; 0.04607; 0.05808];
C0 = [C_ac0; C_isoprop0; C_acetoin0; C_btoh0; C_but0; C_etoh0; C_acetone0];
%Define initial glucose, fructose concentrations
Cgluc(1,1) = C_gluc0;
Cfruc(1,1) = C_fruc0;

%Define fusion parameter
f = 0.00664;  

%% DMMM model with Multi-objective 

for k = 1:480
    %Assign growth efficiency
    if k <104
        bm_eff_clj = 0.3;    
        bm_eff_cac = 0.65;  
    elseif k < 208 && k >= 104              
        bm_eff_clj = 1;
        bm_eff_cac = 0.35;
    else
        bm_eff_clj = 1;
        bm_eff_cac = 0.35;
    end
    %Define total glucose, fructose uptake at each iteration of loop
    ut(1,1) = data(k,2);
    ut(2,1) = data(k,3);
    
    %Exclude non-hybrid Clj from consideration when partitioning glucose
    %consumption between cell type/organisms
    frac(1,1) = [X(k,1)+X(k,3)+X(k,4)+X(k,5)+X(k,6)+X(k,7)+X(k,8)+X(k,9)+X(k,10)+X(k,11)+X(k,2)];
    
    %Partition fructose consumption between all cell type/organisms
    frac(2,1) = sum(X(k,:));
    
    %Initialize summed NGAM by hybrid Cac, Clj, to be used when predicting
    %fermentation yields
    CACfused_ATPuse = 0;
    CLJfused_ATPuse = 0;
    
    %Maximize growth for each organism/cell type included in the simulation
    %(Non-hybrid Cac, Non-hybrid Clj, 5 hybrid metabolism models for both
    %Cac and Clj
    for i = 1:length(species)
        
        %Flux balance analysis for Non-hybrid Cac
        if strcmp(species{i},'Cacpure') == 1
            %load model
            load model.mat
            %Assign substrate utilization based on the fractional abundance
            %of the organism/cell type in the culture
            for j = 1:length(eval(strcat(species{i},'ind')))
                if Cacpureind(j) == 1
                    if X(k,i) == 0
                        modelJoint.lb(4639) = 0; %Cacpure glucose uptake
                        modelJoint.lb(5046) = 0; %Community glucose uptake
                    else
                        modelJoint.lb(4639) = ut(Cacpureind(j))*(X(k,i)/frac(1,1)); %Cacpure glucose uptake
                        modelJoint.lb(5046) = ut(Cacpureind(j))*(X(k,i)/frac(1,1)); %Community glucose uptake
                    end
                else
                    if X(k,i) == 0
                        modelJoint.lb(4557) = 0; %Cacpure fructose uptake
                        modelJoint.lb(5041) = 0; %Community fructose uptake
                    else
                        modelJoint.lb(4557) = ut(Cacpureind(j))*X(k,i); %Cacpure fructose uptake
                        modelJoint.lb(5041) = ut(Cacpureind(j))*X(k,i); %Community fructose uptake
                    end
                end
            end
            
            %Block network activity for all other organisms/cell types in
            %the simulation
            for r = 1476:4550
                modelJoint.lb(r) = 0; %set all other species to zero
                modelJoint.ub(r) = 0; %set all other species to zero
            end
            
            %Assign NGAM based on the fractional abundance the 
            %organism/cell type occupies in the co-culuture
            modelJoint.lb(1100) = 8.39*X(k,i); %Cacpure ATPM
            %Store ATPM for use in fermentation yield FVA
            CACpure_ATPstore = modelJoint.lb(1100);%Cacpure ATPM

            %Perfrom flux balance analysis, store biomass flux
            modelJoint = changeObjective (modelJoint, 'CacpureEx_cpd11416_norm2[u]');
            FBA = optimizeCbModel (modelJoint, 'max');
            mu(k,i) = FBA.f;
            modelJoint.lb(4668) = bm_eff_cac*FBA.f;
            
            %Identify and store feasible ranges for H2, CO2 export by 
            %non-hybrid Cac
            modelJoint = changeObjective (modelJoint, 'CacpureEx_h2[u]');
            h2_store = optimizeCbModel (modelJoint, 'max');
            h2_store_var = h2_store.f;
            modelJoint = changeObjective (modelJoint, 'CacpureEx_co2[u]');
            co2_store = optimizeCbModel (modelJoint, 'max');
            co2_store_var = co2_store.v(4653); %Cac pure CO2 evolution-->Change to [e] ->[u]
            
            %Store non-hybrid Cac glucose and fructose uptake for use in
            %fermentation product FVA
            gluc_flux(k,i) = FBA.v(5046); %Community glucose uptake
            fruc_flux(k,i) = FBA.v(5041); %Community fructose uptake
            
            %X(k,i) = X(k,i)/frac(2,1);
        %Flux balance analysis for Non-hybrid Clj
        elseif strcmp(species{i},'Cljpure') == 1
            load model.mat
            %Assign substrate utilization based on the fractional abundance
            %of the organism/cell type in the culture
            for j = 1:length(eval(strcat(species{i},'ind')))
                if X(k,i) == 0
                    modelJoint.lb(4695) = 0; %Cljpure fructose uptake
                    modelJoint.lb(5041) = 0; %Community fructose uptake
                else
                    modelJoint.lb(4695) = ut(Cljpureind(j))*X(k,i); %Cljpure fructose uptake
                    modelJoint.lb(5041) = ut(Cljpureind(j))*X(k,i); %Community fructose uptake
                end
                modelJoint.lb(4699) = 0; %Cljpure glucose uptake
                modelJoint.lb(5046) = 0; %Community fructose uptake
            end
            %Block non-hybrid Cac flux
            for r = 1:1475
                modelJoint.lb(r) = 0; %set Cac pure bounds to zero
                modelJoint.ub(r) = 0; %set Cac pure bounds to zero
            end
            %Block Cac, Clj hybrid model flux
            for r = 2270:4550
                modelJoint.lb(r) = 0; %set fused model bounds to zero
                modelJoint.ub(r) = 0; %set fused model bounds to zero
            end

            %Place lower bound on Clj H2, CO2 uptake equal to the summed
            %lower bound of H2, CO2 export from all other organsism/cell
            %types in the simulation
            modelJoint.lb(5063) = -h2_store_var; %Cljpure H2 uptake
            modelJoint.lb(4708) = -h2_store_var; %Cljpure H2 uptake
            modelJoint.lb(4686) = -co2_store_var;
            modelJoint.lb(5020) = -co2_store_var;

            %Assign NGAM based on the fractional abundance the 
            %organism/cell type occupies in the co-culuture
            modelJoint.lb(1702) = 0.45*X(k,i); %Cljpure ATPM
            %Store ATPM for use in fermentation yield FVA
            CLJpure_ATPstore = modelJoint.lb(1702); %Cljpure ATPM
            
            %Non-hybrid Clj Flux balance analysis
            modelJoint = changeObjective (modelJoint, 'CljpureBIOMASS_Cl_DSM_WT_46p666M1');
            FBA = optimizeCbModel (modelJoint, 'max');
            
            %Store non-hybrid Clj glucose, fructore uptake, growth rate
            if isempty(FBA.f)
                fruc_flux(k,i) = 0;
                gluc_flux(k,i) = 0;
                mu(k,i) = 0;
            else
                gluc_flux(k,i) = FBA.v(5046); %Community glucose uptake
                fruc_flux(k,i) = FBA.v(5041); %Community fructose uptake
                mu(k,i) = FBA.f;
            end
            
        %Hybrid Cac Flux balance analysis
        elseif strcmp(species{i},'Cacfused') == 1
            load model.mat
            %Set substrate utilization constraints according to relative
            %abundance
            for j = 1:length(eval(strcat(species{i},'ind')))
                if Cacfusedind(j) == 1
                    if XX(k,i) == 0
                        modelJoint.lb(4855) = 0; %Cacfused glucose uptake
                        modelJoint.lb(5046) = 0; %Community glucose uptake
                    else
                        modelJoint.lb(4855) = ut(Cacfusedind(j))*(X(k,i)/frac(1,1)); %Cacfused glucose uptake
                        modelJoint.lb(5046) = ut(Cacfusedind(j))*(X(k,i)/frac(1,1)); %Community glucose uptake
                    end
                else
                     if XX(k,i) == 0
                        modelJoint.lb(4773) = 0; %Cacfused fructose uptake
                        modelJoint.lb(5041) = 0; %Community fructose uptake
                     else
                        modelJoint.lb(4773) = ut(Cacfusedind(j))*X(k,i); %Cacfused fructose uptake
                        modelJoint.lb(5041) = ut(Cacfusedind(j))*X(k,i);  %Community fructose uptake
                     end
                end
            end
            %Block non-hybird model reactions
            for r = 1:2269
                modelJoint.lb(r) = 0; %set non-hybrid models bounds to zero
                modelJoint.ub(r) = 0; %set non-hybrid models bounds to zero
            end
            %Block hybrid Clj model reactions
            for r = 3752:4550
                modelJoint.lb(r) = 0; %set hybrid Clj bounds to zero
                modelJoint.ub(r) = 0; %set hybrid Clj bounds to zero
            end
            %Assign NGAM based on the fractional abundance the 
            %organism/cell type occupies in the co-culuture
            modelJoint.lb(3369) = 8.39*X(k,i); %Cac hybrid ATPM
            %Store ATPM for use in fermentation yield FVA
            CACfused_ATPstore = modelJoint.lb(3369); %Cac hybrid ATPM
            CACfused_ATPuse = CACfused_ATPuse + CACfused_ATPstore; 
            %Hybrid Cac Flux balance analysis
            modelJoint = changeObjective (modelJoint, 'CacfusedEx_cpd11416_norm2[u]');
            FBA = optimizeCbModel (modelJoint, 'max');
            %Store glucose, fructose uptake and growth for use during 
            %Euler integration, fermentation product FVA
            gluc_flux(k,i) = FBA.v(5046); %Community glucose uptake
            fruc_flux(k,i) = FBA.v(5041); %Community fructose uptake
            mu(k,i) = FBA.f;
        
        %Hybrid Clj Flux balance analysis   
        elseif strcmp(species{i},'Cljfused') == 1
            load model.mat
            %Set substrate utilization constraints according to relative
            %abundance
            for j = 1:length(eval(strcat(species{i},'ind')))
                if Cljfusedind(j) == 1
                    if XX(k,i) == 0
                        modelJoint.lb(4917) = 0 %Clj fused glucose uptake
                        modelJoint.lb(5046) = 0 %Community glucose uptake
                    else
                        modelJoint.lb(4917) = ut(Cljfusedind(j))*X(k,i); %Clj fused glucose uptake
                        modelJoint.lb(5046) = ut(Cljfusedind(j))*X(k,i); %Community glucose uptake
                    end
                else
                    if XX(k,i) == 0
                        modelJoint.lb(4913) = 0; %Clj fused fructose uptake
                        modelJoint.lb(5041) = 0; %Community fructose uptake
                    else    
                        modelJoint.lb(4913) = ut(Cljfusedind(j))*X(k,i); %Clj fused fructose uptake
                        modelJoint.lb(5041) = ut(Cljfusedind(j))*X(k,i); %Community fructose uptake
                    end
                end
            end
            %Block non-hybird model and hybrid Cac reactions
            for r = 1:3751
                modelJoint.lb(r) = 0; %set bounds for all other models to zero
                modelJoint.ub(r) = 0; %set bounds for all other models to zero
            end
            %Assign NGAM based on the fractional abundance the 
            %organism/cell type occupies in the co-culuture
            modelJoint.lb(3978) = 0.45*X(k,i); %Cljfused ATPM
            %Store ATPM for use in fermentation yield FVA
            CLJfused_ATPstore = modelJoint.lb(3978); %Cljfused ATPM
            CLJfused_ATPuse = CLJfused_ATPuse + CLJfused_ATPstore;
            
            %Hybrid Clj flux balance analysis
            modelJoint = changeObjective (modelJoint, 'CljfusedBIOMASS_Cl_DSM_WT_46p666M1');
            FBA = optimizeCbModel (modelJoint, 'max');
            %Store glucose, fructose uptake and growth for use during 
            %Euler integration, fermentation product FVA
            gluc_flux(k,i) = FBA.v(5046);
            fruc_flux(k,i) = FBA.v(5041);
            mu(k,i) = FBA.f;
        end

    end
    
    %Determine summed growth, H2 and CO2 cross-feeding for Cac hybrid
    %states
    if (mu(k,2) + mu(k,4) + mu(k,6) + mu(k,8) + mu(k,10))-(mu(k,2) + mu(k,4) + mu(k,6) + mu(k,8) + mu(k,10))*0.01 > 0
        load model.mat
            
        for r = 1:2269
           modelJoint.lb(r) = 0; %set non-hybrid models bounds to zero
           modelJoint.ub(r) = 0; %set non-hybrid models bounds to zero
        end
        for r = 3752:4550
           modelJoint.lb(r) = 0; %set hybrid Clj bounds to zero
           modelJoint.ub(r) = 0; %set hybrid Clj bounds to zero
        end
        modelJoint.lb(5046) = gluc_flux(k,2) + gluc_flux(k,4) + gluc_flux(k,6) + gluc_flux(k,8) + gluc_flux(k,10);
        modelJoint.lb(5041) = fruc_flux(k,2) + fruc_flux(k,4) + fruc_flux(k,6) + fruc_flux(k,8) + fruc_flux(k,10);
        modelJoint.lb(4773) = fruc_flux(k,2) + fruc_flux(k,4) + fruc_flux(k,6) + fruc_flux(k,8) + fruc_flux(k,10);
        modelJoint.lb(4855) = gluc_flux(k,2) + gluc_flux(k,4) + gluc_flux(k,6) + gluc_flux(k,8) + gluc_flux(k,10);
        
        modelJoint.lb(3369) = 8.39*(X(k,2) + X(k,4) + X(k,6) + X(k,8) + X(k,10));
        
        modelJoint = changeObjective (modelJoint, 'CacfusedEx_cpd11416_norm2[u]');
        FBA = optimizeCbModel (modelJoint, 'max');
        modelJoint.lb(4884) = bm_eff_cac*FBA.f;
        store_fuse_mu = bm_eff_cac*FBA.f;
    
        modelJoint = changeObjective (modelJoint, 'ExCom_h2[u]');
        FBA = optimizeCbModel (modelJoint, 'max');
        h2_plug = FBA.f;
        h2_store_var = h2_plug + h2_store_var;
        modelJoint = changeObjective (modelJoint, 'ExCom_co2[u]');
        FBA = optimizeCbModel (modelJoint, 'max');
        co2_plug = FBA.f;
        co2_store_var = co2_plug + co2_store_var;
        % save 
        co2_Cac{k,1} = co2_store_var;
        h2_Cac{k,1} = h2_store_var;
    end
    
    %Determine summed growth, H2 and CO2 cross-feeding for Clj hybrid
    %states
    if (mu(k,3) + mu(k,5) + mu(k,7) + mu(k,9) + mu(k,11))-(mu(k,3) + mu(k,5) + mu(k,7) + mu(k,9) + mu(k,11))*0.01 > 0
        load model.mat
        for r = 1:3751
                modelJoint.lb(r) = 0; %set bounds for all other models to zero
                modelJoint.ub(r) = 0; %set bounds for all other models to zero
        end
        
        modelJoint.lb(5046) = gluc_flux(k,3) + gluc_flux(k,5) + gluc_flux(k,7) + gluc_flux(k,9) + gluc_flux(k,11);
        modelJoint.lb(5041) = fruc_flux(k,3) + fruc_flux(k,5) + fruc_flux(k,7) + fruc_flux(k,9) + fruc_flux(k,11);
        modelJoint.lb(4913) = fruc_flux(k,3) + fruc_flux(k,5) + fruc_flux(k,7) + fruc_flux(k,9) + fruc_flux(k,11);
        modelJoint.lb(4917) = gluc_flux(k,3) + gluc_flux(k,5) + gluc_flux(k,7) + gluc_flux(k,9) + gluc_flux(k,11);
        
        modelJoint.lb(3978) = 0.45*(X(k,3) + X(k,5) + X(k,7) + X(k,9) + X(k,11));
        
        modelJoint = changeObjective (modelJoint, 'CljfusedBIOMASS_Cl_DSM_WT_46p666M1');
        FBA = optimizeCbModel (modelJoint, 'max');
        modelJoint.lb(3752) = bm_eff_clj*FBA.f;
        store_fuse_Clj_mu = bm_eff_clj*FBA.f;
    
        modelJoint = changeObjective (modelJoint, 'ExCom_h2[u]');
        FBA = optimizeCbModel (modelJoint, 'max');
        h2_plug_Clj = FBA.f;
        h2_store_var = h2_plug_Clj + h2_store_var;
        modelJoint = changeObjective (modelJoint, 'ExCom_co2[u]');
        FBA = optimizeCbModel (modelJoint, 'max');
        co2_plug_Clj = FBA.f;
        if ~isempty(co2_plug_Clj)
            co2_store_var = co2_plug_Clj + co2_store_var;
        else
            co2_store_var = co2_store_var;
        end
            %Store Clj CO2, H2 uptake rates
    co2_Clj{k,1} = co2_store_var;
    h2_Clj{k,1} = h2_store_var;
    
    end
    
    load model.mat    
    
    %Set lower bounds for community glucose and fructose uptake rates
    modelJoint.lb(5046) = sum(gluc_flux(k,:));
    modelJoint.lb(5041) = sum(fruc_flux(k,:));
    
    %Set lower bounds for organism/cell type glucose, fructose uptake rates
    modelJoint.lb(4639) = gluc_flux(k,1);
    modelJoint.lb(4699) = gluc_flux(k,12);
    modelJoint.lb(4855) = gluc_flux(k,2) + gluc_flux(k,4) + gluc_flux(k,6) + gluc_flux(k,8) + gluc_flux(k,10);
    modelJoint.lb(4917) = gluc_flux(k,3) + gluc_flux(k,5) + gluc_flux(k,7) + gluc_flux(k,9) + gluc_flux(k,11);
    modelJoint.lb(4557) = fruc_flux(k,1); 
    modelJoint.lb(4695) = fruc_flux(k,12);
    modelJoint.lb(4773) = fruc_flux(k,2) + fruc_flux(k,4) + fruc_flux(k,6) + fruc_flux(k,8) + fruc_flux(k,10);
    modelJoint.lb(4913) = fruc_flux(k,3) + fruc_flux(k,5) + fruc_flux(k,7) + fruc_flux(k,9) + fruc_flux(k,11);
    
    %Place lower bound on organism/cell type growth rate (99% of maximum
    %growth rate)
    modelJoint.lb(4668) = bm_eff_cac*mu(k,1)-bm_eff_cac*mu(k,1)*0.01;
    modelJoint.lb(1476) = bm_eff_clj*mu(k,12)-bm_eff_clj*mu(k,12)*0.01;
    if (mu(k,2) + mu(k,4) + mu(k,6) + mu(k,8) + mu(k,10))-(mu(k,2) + mu(k,4) + mu(k,6) + mu(k,8) + mu(k,10))*0.01 > 0
        modelJoint.lb(4884) = store_fuse_mu - 0.01*store_fuse_mu;
    else
        0;
    end
    if (mu(k,3) + mu(k,5) + mu(k,7) + mu(k,9) + mu(k,11))-(mu(k,3) + mu(k,5) + mu(k,7) + mu(k,9) + mu(k,11))*0.01 > 0
        modelJoint.lb(3752) = store_fuse_Clj_mu - 0.01*store_fuse_Clj_mu;
    else
        0;
    end

% objective sets for the single objective case
%     %Define fermentation product export flux as FBA objective
%     if k <= 119 
%         rxnsList = {'ExCom_but[u]'}; 
%     elseif k > 119 && k<313
%         rxnsList = {'ExCom_lac__L[u]'};  
%     else
%         rxnsList = {'ExCom_lac__L[u]'};
%     end
    modelJoint = changeObjective (modelJoint, rxnsList{1}); 

    %Block flux through hybrid states with zero relative abundance
    if X(k,1) == 0
        for r = 1:1475
            modelJoint.lb(r) = 0;
            modelJoint.ub(r) = 0;
        end
    end
    if X(k,12) == 0
        for r = 1476:2269
            modelJoint.lb(r) = 0;
            modelJoint.ub(r) = 0;
        end
    end
    if X(k,2) + X(k,4) + X(k,6) + X(k,8) + X(k,10) == 0
        for r = 2270:3751
            modelJoint.lb(r) = 0;
            modelJoint.ub(r) = 0;
        end
    end
    if X(k,3) + X(k,5) + X(k,7) + X(k,9) + X(k,11) == 0
        for r = 3746:4550
            modelJoint.lb(r) = 0;
            modelJoint.ub(r) = 0;
        end
    end
        
    %Place lower bounds on Clj H2 and CO2 uptake rate
    modelJoint.lb(4708) = -h2_store_var;
    modelJoint.lb(4686) = -co2_store_var; 
    %Store Clj CO2, H2 uptake rates
    co2_vals{k,1} = co2_store_var;
    h2_vals{k,1} = h2_store_var;
    
    %Place NGAM lower bounds for each organism/cell type
    modelJoint.lb(1100) = CACpure_ATPstore;
    modelJoint.lb(1702) = CLJpure_ATPstore;
    modelJoint.lb(3369) = CACfused_ATPuse;
    modelJoint.lb(3978) = CLJfused_ATPuse;
    
    %Place bounds on exchange fluxes to preciesly match flux balance
    %analysis results with experimental fermentation profiles
    modelJoint.lb(4985) = data(k,5); %isopropanol
    modelJoint.lb(4987) = data(k,4); %acetate
    modelJoint.lb(4988) = data(k,10); %acetone
    modelJoint.lb(4991) = data(k,6); %acetoin
    modelJoint.lb(5034) = data(k,9); %ethanol
    modelJoint.lb(5011) = data(k,7); %butanol

% multi-objective function derived from the integrated approach
    if k < 80 
       f=[zeros(4984,1);0.14;0;0.3; 0 ;zeros(22,1);0.112;0.287;zeros(21,1);0.161;zeros(107,1)]';       
       modelJoint.lb(5012) =  0;  

      elseif k >= 80  && k< 110  
         f=[zeros(4984,1);0.078; 0 ;0.13; 0 ; zeros(22,1);0.128 ; 0.187 ;zeros(21,1);0.477;zeros(107,1)]';
    modelJoint.lb(5012) = data(k,8);

    elseif k >= 110  && k< 130   
         f=[zeros(4984,1);0.047; 0 ;0.169; 0 ;zeros(22,1);0.073;0.179;zeros(21,1);0.532;zeros(107,1)]';

    elseif k >= 130  && k< 170   
         f=[zeros(4984,1);0.22; 0 ;0.164; 0.166 ; zeros(22,1);0.281; 0 ;zeros(21,1);0.169;zeros(107,1)]';

    elseif k >= 170 && k<210  
        f=[zeros(4984,1);0.228; 0 ;0.192;0.151;zeros(22,1);0.192;0;zeros(21,1);0.237;zeros(107,1)]';

    elseif k >= 210  && k<270 
        f=[zeros(4984,1);0.181;0;0.054;0;zeros(22,1);0.216;0;zeros(21,1);0.549;zeros(107,1)]';
    
    elseif k >= 270 && k<310  
        f=[zeros(4984,1);0.588;0;0.03;0;zeros(22,1);0.268;0;zeros(21,1);0.114;zeros(107,1)]';

    elseif k >= 310 && k<323 
        f=[zeros(4984,1);0.216;0;0.159;0;zeros(22,1);0.243;0;zeros(21,1);0.382;zeros(107,1)]';

    else k >= 323 
        f=[zeros(4984,1);0.216;0;0.159;0;zeros(22,1);0.243;0;zeros(21,1);0.382;zeros(107,1)]';

    end

    %Estimate fermentation yields
    % inequality
        struc=full(modelJoint.S);
        LB=modelJoint.lb;
        UB=modelJoint.ub;
        b=(modelJoint.b);
     
    FBA_multi = linprog(f,[],[],struc,b,LB,UB);

%    FBA = optimizeCbModel (modelJoint, 'max');

    save_lb=modelJoint.lb;
    SLB=[SLB,save_lb];
    save_ub=modelJoint.ub;
    SUB=[SUB,save_ub];

    %Store solution
    if ~isempty( FBA_multi )

        store_obj(k,1) =  FBA_multi (4987); %acetate
        store_obj(k,2) =  FBA_multi (4985); %isopropanol
        store_obj(k,4) =  FBA_multi (5011); %butanol
        store_obj(k,6) =  FBA_multi (5034); %ethanol

        if k == 1
            for rr = 1:4
                dStore_obj(k,rr) = store_obj(k,rr)*dt;
                Cstore_obj(k,rr) = C0(rr) + dStore_obj(k,rr);
                Cstore_GL(k,rr) = Cstore_obj(k,rr) * MW(rr); %concentration g/L
            end
        else
            for rr = 1:4
                dStore_obj(k,rr) = store_obj(k,rr)*dt;
                Cstore_obj(k,rr) = Cstore_obj(k-1,rr) + dStore_obj(k,rr);
                Cstore_GL(k,rr) = Cstore_obj(k,rr) * MW(rr); %concentration g/L
            end
        end
    end
end


%% plot the predicted results

subplot(2,2,1)
plot(concentration(1:5,1),concentration(1:5,7)*0.0601,'*')
hold on
xs=linspace(0,0.1*(k),k);
plot(xs,Cstore_GL(1:k,2))
title('isopropanol')
xlabel('Time (hr)')
ylabel('Titer (g/L)')
subplot(2,2,2)
plot(concentration(1:5,1),concentration(1:5,10)*0.074121,'*')
hold on
xs=linspace(0,0.1*(k),k);
plot(xs,Cstore_GL(1:k,4))
title('butanol')
xlabel('Time (hr)')
ylabel('Titer (g/L)')
subplot(2,2,3)
plot(concentration(1:5,1),concentration(1:5,6)*0.04607,'*')
hold on
xs=linspace(0,0.1*(k),k);
plot(xs,Cstore_GL(1:k,6))
title('ethanol')
xlabel('Time (hr)')
ylabel('Titer (g/L)')
subplot(2,2,4)
plot(concentration(1:5,1),concentration(1:5,4)*0.060052,'*')
hold on
xs=linspace(0,0.1*(k),k);
plot(xs,Cstore_GL(1:k,1))
title('acetate')
xlabel('Time (hr)')
ylabel('Titer (g/L)')

