% ref: Au, J., Choi, J., Jones, S. W., Venkataramanan, K. P., & Antoniewicz, M. R. (2014). Parallel labeling experiments validate Clostridium acetobutylicum metabolic network model for 13C metabolic flux analysis. Metabolic engineering, 26, 23-33.
% model = buildModelForMFG_fromMFA(excelPath); % excel path to the data file
% or
load data_example

% 2) Use the MFA best-fit vector
v = model.vMFA;

% 3) Build the MFG
mfg = buildMFG_fromFlux(model, v);

% 4) Pick source/sink reactions by (partial) name
rxnNames = string(model.rxnNames);

% Source: hexokinase step (~100), exactly the MFA normalization
jSrc = find(contains(rxnNames,'Gluc.Ext + ATP -> G6P + ADP'),1);

% Sinks: product secretions (adjust if your strings differ)
jBut = find(contains(rxnNames,'Bt -> Bt.Ext'),1);      % butyrate/butanol? (use name in your list)
jAce = find(contains(rxnNames,'Acne -> Acne.Ext'),1);  % acetone
jEt  = find(contains(rxnNames,'EtOH -> EtOH.Ext'),1);  % ethanol
jAc  = find(contains(rxnNames,'Ac -> Ac.Ext'),1);      % acetate

% 5) Max-flow or path visualization
G = mfg.G;
figure; H = plot(G, 'Layout','force', 'MarkerSize',3); title('MFG');
highlight(H, jSrc, 'NodeColor','b','MarkerSize',8);

if ~isempty(jBut)
    [mf,GF,~,~] = maxflow(G, jSrc, jBut);
    highlight(H, GF, 'EdgeColor','g','LineWidth',2);
    highlight(H, jBut, 'NodeColor','g','MarkerSize',8);
    fprintf('Max-flow Glucose→But product = %.3f (arbitrary MFG units)\n', mf);
end


% Built model: 92 mets x 81 rxns; vMFA length = 81 (aligned)
% Max-flow Glucose→But product = 64.601 (arbitrary MFG units)

%%
S  = full(model.S);
r  = S*model.vMFA - model.b;              % residual by metabolite
absr = abs(r);
[absr_sorted, idx] = sort(absr, 'descend');

disp('Top residual metabolites:');
for k = 1:min(15,numel(idx))
    fprintf('%2d) %-15s  residual = %+8.3f\n', k, model.mets{idx(k)}, r(idx(k)));
end


%%

rxnNames = string(model.rxnNames);
G = mfg.G;

% example source/sink
jSrc = find(contains(rxnNames,'Gluc.Ext + ATP -> G6P + ADP'),1);
jBut = find(contains(rxnNames,'Bt -> Bt.Ext'),1);

[mf,GF] = maxflow(G, jSrc, jBut);
fprintf('Max-flow = %.3f\n', mf);

% map end-node names to numeric indices
sIdx = findnode(G, GF.Edges.EndNodes(:,1));
tIdx = findnode(G, GF.Edges.EndNodes(:,2));
w    = GF.Edges.Weight;

tbl = table(rxnNames(sIdx), rxnNames(tIdx), w, ...
    'VariableNames', {'from','to','max-flow'});
disp(tbl)

%%

rxnNames = string(model.rxnNames);
G = mfg.G;

% Source = hexokinase step (your scale anchor)
jSrc = find(contains(rxnNames,'Gluc.Ext + ATP -> G6P + ADP'),1);

% Define sinks and the matching EX patterns in your rxnNames
sinkList = {
  'butanol',  'BtOH -> BtOH.Ext';
  'butyrate', 'Bt -> Bt.Ext';
  'acetone',  'Acne -> Acne.Ext';
  'ethanol',  'EtOH -> EtOH.Ext';
  'acetate',  'Ac -> Ac.Ext';
};

nS = size(sinkList,1);
mfVals  = nan(nS,1);    % MFG maxflow
mfaVals = nan(nS,1);    % MFA EX flux (Best Fit)
names   = strings(nS,1);

for i = 1:nS
    names(i) = sinkList{i,1};
    pat      = sinkList{i,2};

    jSink = find(contains(rxnNames, pat), 1);
    if isempty(jSink)
        warning('Sink not found: %s', pat);
        continue
    end

    % MFG max-flow Glucose -> sink
    [mf,~] = maxflow(G, jSrc, jSink);
    mfVals(i) = mf;

    % MFA value for the same EX reaction
    mfaVals(i) = model.vMFA(jSink);
end

% Summary table
Tcmp = table(names, mfVals, mfaVals, mfaVals./model.vMFA(jSrc), ...
    'VariableNames', {'sink','MFG','MFA','MFA_over_Gluc'});
disp(Tcmp)

% Correlations (guard n>=2)
ok = isfinite(mfVals) & isfinite(mfaVals);
n  = nnz(ok);
if n >= 2
    rhoS = corr(mfVals(ok), mfaVals(ok), 'type','Spearman');
    R2   = corr(mfVals(ok), mfaVals(ok))^2;
    fprintf('MFG vs MFA on sinks: Spearman=%.3f, R^2=%.3f (n=%d)\n', rhoS, R2, n);
else
    fprintf('Not enough sinks to compute correlation (n=%d)\n', n);
end


%%
% Example: a simple objective – maximize combined solvents (replace with your TIOBJFIND c*)
cStar = zeros(size(model.S,2),1);
for pat = ["BtOH -> BtOH.Ext","Bt -> Bt.Ext","Acne -> Acne.Ext","EtOH -> EtOH.Ext","Ac -> Ac.Ext"]
    j = find(contains(rxnNames, pat),1);
    if ~isempty(j), cStar(j) = 1; end
end

% Solve LP: maximize c*'v
f = -cStar;  % linprog minimizes
opts = optimoptions('linprog','Display','none');
[v_FBA, fval, exitflag] = linprog(f,[],[],full(model.S),model.b,model.lb,model.ub,opts);
assert(exitflag>0,'FBA LP infeasible with chosen bounds/objective — check HK & EX bounds');

% 1) Make feasible
model = makeFeasibleForObjective(model);

% 2) Anchor the HK step (scale) ≈ 100
rxnNames = string(model.rxnNames);
jHK = find(contains(rxnNames,'Gluc.Ext + ATP -> G6P + ADP'),1);
assert(~isempty(jHK),'HK step not found in model.rxnNames');
model.lb(jHK) =  99.9;
model.ub(jHK) = 100.1;

% 3) Define a simple objective (example: sum of product secretions)
cStar = zeros(size(model.S,2),1);
for pat = ["Bt -> Bt.Ext","BtOH -> BtOH.Ext","Acne -> Acne.Ext","EtOH -> EtOH.Ext","Ac -> Ac.Ext"]
    j = find(contains(rxnNames, pat),1);
    if ~isempty(j), cStar(j) = 1; end
end
if ~any(cStar)
    warning('No product exchanges matched; using HK as placeholder objective');
    cStar(jHK) = 1;
end

% 4) Solve LP: maximize cStar''*v  (linprog minimizes, so use -cStar)
f = -cStar;
opts = optimoptions('linprog','Display','final');  % show status
[Sfull, b] = deal(full(model.S), model.b);
[v_FBA, fval, exitflag] = linprog(f,[],[],Sfull,b,model.lb,model.ub,opts);

fprintf('linprog exitflag = %d\n', exitflag);
assert(exitflag>0,'FBA LP infeasible with chosen bounds/objective — check HK & EX bounds');

% 5) Report key numbers
fprintf('Objective value (max) = %.4f\n', -fval);
resid = norm(Sfull*v_FBA - b);
fprintf('Mass-balance residual ||S*v - b|| = %.2e\n', resid);

% Show a few exchange fluxes
show = ["Gluc.Ext + ATP -> G6P + ADP", ...
        "Ac -> Ac.Ext","Bt -> Bt.Ext","Acne -> Acne.Ext","EtOH -> EtOH.Ext","BtOH -> BtOH.Ext","CO2 -> CO2.Ext","H2 -> H2.Ext"];
for p = show
    j = find(contains(rxnNames,p),1);
    if ~isempty(j)
        fprintf('v(%s) = %.4f  [%.1f, %.1f]\n', p, v_FBA(j), model.lb(j), model.ub(j));
    end
end

%%
%% replace with CoI c*
% ---- map reaction names -> indices in your model ----
rxnNames = string(model.rxnNames);

jBtOH = find(contains(rxnNames,'BtOH -> BtOH.Ext'),1);   % butanol
jBt   = find(contains(rxnNames,'Bt -> Bt.Ext'),1);       % butyrate
jAcne = find(contains(rxnNames,'Acne -> Acne.Ext'),1);   % acetone
jEtOH = find(contains(rxnNames,'EtOH -> EtOH.Ext'),1);   % ethanol
jAc   = find(contains(rxnNames,'Ac -> Ac.Ext'),1);       % acetate

% jbiomass = find(contains(rxnNames,'0.410 Ala + 0.070 Arg + 0.082 Asn + 0.082 Asp +  … -> Biomass +'),1);   % butanol


assert(~isempty(jBtOH) && ~isempty(jBt) && ~isempty(jAcne) ...
    && ~isempty(jEtOH) && ~isempty(jAc), 'One or more EX sinks not found.');

% Build split model (forward/backward) once:
S = full(model.S); [m,n] = size(S);
rev = model.lb < 0;       % reversible flags
S2  = [ S, -S(:,rev) ];   % v+ for all; v- only for reversible
lb2 = [ max(model.lb,0) ; zeros(nnz(rev),1) ];   % v+ >= 0; v- >= 0
ub2 = [ model.ub         ; model.ub(rev)   ];
b   = model.b;


% 13C BestFit (your MFA numbers)
mfa.butOH = 5.0718;   % BtOH
mfa.but   = 64.601;   % Bt
mfa.acetone = 2.2874; % Acetone
mfa.etoh  = 1.829;    % EtOH
mfa.ac    = 30.575;   % Acetate

% find indices (already done earlier)
rxnNames = string(model.rxnNames);
jBtOH = find(contains(rxnNames,'BtOH -> BtOH.Ext'),1);
jBt   = find(contains(rxnNames,'Bt -> Bt.Ext'),1);
jAcne = find(contains(rxnNames,'Acne -> Acne.Ext'),1);
jEtOH = find(contains(rxnNames,'EtOH -> EtOH.Ext'),1);
jAc   = find(contains(rxnNames,'Ac -> Ac.Ext'),1);

% (1) Anchor HK ~100
jHK = find(contains(rxnNames,'Gluc.Ext + ATP -> G6P + ADP'),1);
model.lb(jHK)=99.9; model.ub(jHK)=100.1;

% (2) Set product bounds to 13C BestFit ± small tolerance (e.g., 5%)
tol = 0.05;
model.lb(jBtOH) = max(0, (1-tol)*mfa.butOH);  model.ub(jBtOH) = (1+tol)*mfa.butOH;
model.lb(jBt)   = max(0, (1-tol)*mfa.but);    model.ub(jBt)   = (1+tol)*mfa.but;
model.lb(jAcne) = max(0, (1-tol)*mfa.acetone);model.ub(jAcne) = (1+tol)*mfa.acetone;
model.lb(jEtOH) = max(0, (1-tol)*mfa.etoh);   model.ub(jEtOH) = (1+tol)*mfa.etoh;
model.lb(jAc)   = max(0, (1-tol)*mfa.ac);     model.ub(jAc)   = (1+tol)*mfa.ac;

% (3) Build c from your weights (any normalization scheme is fine)
cStar = zeros(size(model.S,2),1);
% e.g., sum-normalized weights from your MFG:
cStar(jBtOH)=0.0485971243; cStar(jBt)=0.6189957859; cStar(jAcne)=0.0219174774;
cStar(jEtOH)=0.0175251667; cStar(jAc)=0.2929644457;

% (4) Solve FBA
f = -cStar;
opts = optimoptions('linprog','Display','none');
[v_FBA,fval,exitflag] = linprog(f,[],[],full(model.S),model.b,model.lb,model.ub,opts);
assert(exitflag>0,'LP infeasible; relax tolerances slightly.');

% (5) Print achieved product fluxes vs 13C targets
fprintf('BtOH %.4f (target %.4f)\n', v_FBA(jBtOH), mfa.butOH);
fprintf('Bt   %.4f (target %.4f)\n', v_FBA(jBt),   mfa.but);
fprintf('Acetone %.4f (target %.4f)\n', v_FBA(jAcne), mfa.acetone);
fprintf('EtOH %.4f (target %.4f)\n',   v_FBA(jEtOH), mfa.etoh);
fprintf('Acetate %.4f (target %.4f)\n',v_FBA(jAc),   mfa.ac);
%%


mfa = struct('butOH',5.0718,'but',64.6010,'acetone',2.2874,'etoh',1.8290,'ac',30.5750);
vals = [v_FBA(jBtOH) v_FBA(jBt) v_FBA(jAcne) v_FBA(jEtOH) v_FBA(jAc)];
targ = [mfa.butOH   mfa.but    mfa.acetone   mfa.etoh     mfa.ac];
names = {'BtOH','Bt','Acetone','EtOH','Acetate'};
for i=1:numel(names)
    fprintf('%-7s: %.4f (target %.4f) ', ...
        names{i}, vals(i), targ(i), 100*(vals(i)-targ(i))/targ(i));
end

%% 
% using TIObjFind CoI
% BtOH 5.3254 (target 5.0718)
% Bt   67.8311 (target 64.6010)
% Acetone 2.4018 (target 2.2874)
% EtOH 1.9204 (target 1.8290)
% Acetate 32.1037 (target 30.5750)

% using biomass rxn as the objective
% BtOH   : 4.8182 (target 5.0718)
% Bt     : 61.3709 (target 64.6010) 
% Acetone: 2.1730 (target 2.2874) 
% EtOH   : 1.9204 (target 1.8290) 
% Acetate: 29.0462 (target 30.5750)


