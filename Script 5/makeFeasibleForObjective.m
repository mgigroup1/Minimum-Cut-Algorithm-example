function model = makeFeasibleForObjective(model)
% Adds exchange rxns for every [e] metabolite if missing.
% Initializes/extends model.c only if it already exists.

    rxnNames = string(model.rxnNames);
    mets     = string(model.mets);
    nMets0   = numel(mets);
    nRxn0    = size(model.S,2);

    % 0) add EX for every [e] metabolite not already present
    addCols = {};
    addNames = {};
    addIDs   = {};
    for i = 1:nMets0
        mi = mets(i);
        if ~endsWith(mi,'[e]'), continue; end
        exName = "Exchange " + mi;                 % name
        if any(contains(rxnNames, exName)), continue; end

        exCol = sparse(nMets0,1);  exCol(i) = -1;  % X[e] -> ∅
        addCols{end+1}  = exCol;                   %#ok<AGROW>
        addNames{end+1} = char(exName);
        base = erase(mi,'[e]');
        addIDs{end+1}   = char("EX_"+base+"(e)");
    end

    if ~isempty(addCols)
        model.S        = [model.S, cat(2,addCols{:})];
        model.lb       = [model.lb;  -1000*ones(numel(addCols),1)];
        model.ub       = [model.ub;  +1000*ones(numel(addCols),1)];
        if isfield(model,'b'); model.b = model.b; else, model.b = zeros(size(model.S,1),1); end
        model.rxns     = [model.rxns;     addIDs(:)];
        model.rxnNames = [model.rxnNames; addNames(:)];
        % extend c only if it exists
        if isfield(model,'c')
            model.c = [model.c; zeros(numel(addCols),1)];
        end
        % extend vMFA if you keep it in the same struct
        if isfield(model,'vMFA')
            model.vMFA = [model.vMFA; NaN(numel(addCols),1)];
        end
    end

    % 1) Allow uptake/secretion where reasonable
    rxnNames = string(model.rxnNames);  % refresh
    allowUptake   = ["Gluc[e]","NH3[e]","SO4[e]"];
    allowSecretion= ["CO2[e]","H2[e]","Ac[e]","Bt[e]","EtOH[e]","Acne[e]","BtOH[e]"];

    for k = 1:numel(allowUptake)
        j = find(contains(rxnNames,"Exchange "+allowUptake(k)),1);
        if ~isempty(j), model.lb(j) = -1000; model.ub(j) = +1000; end
    end
    for k = 1:numel(allowSecretion)
        j = find(contains(rxnNames,"Exchange "+allowSecretion(k)),1);
        if ~isempty(j), model.lb(j) = -1000; model.ub(j) = +1000; end
    end
end
