function model = buildModelFromMFA(excelPath)
% Build a COBRA-style model struct from Au et al. (mmc1.xlsx)
% model.S (mets x rxns), model.rxns, model.rxnNames, model.mets,
% model.lb, model.ub, model.c, model.b, model.description
%
% Notes:
% - Uses the Minimal Model columns (BestFit/LB95/UB95).
% - Directionality: '->' treated irreversible (lb = 0 by default),
%                   '<=>' treated reversible (lb = -1000 by default).
% - Bounds: by default sets wide bounds; optionally narrows around MFA flux.
%   See the "TUNE BOUNDS" section to choose your preference.

    if nargin < 1
        error('Provide path to 1-s2.0-S1096717614001050-mmc1.xlsx');
    end

    % -------- 1) Read and locate the main table  --------
    raw = readcell(excelPath, 'Sheet', 'Results of 13C-MFA');

    % Find the header row that starts the intracellular fluxes table
    hdrRow = find(strcmp(raw(:,1),'INTRACELLULAR FLUXES'), 1, 'first');
    assert(~isempty(hdrRow), 'Could not find "INTRACELLULAR FLUXES" header.');
    % The actual column headers are a couple rows below
    colHdrRow = hdrRow + 1;

    % Expected header layout after the header row:
    % Cols: [Flux No., Reaction, Ext_Best, Ext_LB95, Ext_UB95, <blank>, Min_Best, Min_LB95, Min_UB95]
    % We’ll just pull A:I and rename columns.
    T = cell2table(raw(colHdrRow+1:end,1:9), 'VariableNames', ...
        {'FluxNo','Reaction','ExtBest','ExtLB','ExtUB','Blank','MinBest','MinLB','MinUB'});

    % Drop all-NaN rows and rows before data begins
    % Keep rows where FluxNo is numeric OR Reaction is text
    isData = cellfun(@(x) isnumeric(x) && ~isnan(x), T.FluxNo) | ...
             cellfun(@(x) ischar(x) || isstring(x), T.Reaction);
    T = T(isData,:);

    % Many header/blank lines still sneak in; keep only rows where Reaction looks like a reaction
    looksLikeRxn = cellfun(@(s) ~isempty(s) && contains(string(s), {'->','<=>'}), T.Reaction);
    T = T(looksLikeRxn,:);

    % Minimal model presence: if the row in raw had "not present in model" in the minimal columns,
    %_minimalBest will be NaN or the Reaction row above said so. Just filter on MinBest numeric.
    keepMinimal = cellfun(@(x) isnumeric(x) && ~isnan(x), T.MinBest);
    T = T(keepMinimal,:);

    % -------- 2) Parse reactions -> stoichiometry  --------
    rxnNames = string(T.Reaction);
    nRxn = numel(rxnNames);
    mets = string([]);                 % grow list of metabolite IDs
    metToIdx = containers.Map('KeyType','char','ValueType','int32');

    % Helper to add/get metabolite index
    function idx = addMet(met)
        met = strtrim(met);
        if ~isKey(metToIdx, met)
            metToIdx(met) = length(mets) + 1;
            mets(end+1) = met; %#ok<AGROW>
        end
        idx = metToIdx(met);
    end

    % Build triplets for sparse S
    I = []; J = []; V = [];

    arrows = strings(nRxn,1);
    rxnIDs = strings(nRxn,1);

    for j = 1:nRxn
        rstr = rxnNames(j);

        if contains(rstr,'<=>')
            arrows(j) = '<=>';
            parts = split(rstr,'<=>');
        elseif contains(rstr,'->')
            arrows(j) = '->';
            parts = split(rstr,'->');
        else
            error('Reaction %d has no arrow.', j);
        end
        if numel(parts) ~= 2
            error('Reaction %d parse error: %s', j, rstr);
        end

        left = strtrim(parts(1));
        right = strtrim(parts(2));

        % Strip trailing annotations like "(net)" or "(exch)" on either side
        left  = erase(left,  {'(net)','(exch)'});
        right = erase(right, {'(net)','(exch)'});

        % Tokenize each side by '+'
        leftTerms  = split(left,  '+');
        rightTerms = split(right, '+');

        % Parse a term of the form "[coef] MetName"
        % coef optional, may be decimal; MetName alnum/._- (e.g., Gluc.Ext, E-C2)
        parseTerm = @(s) localParseTerm(s);

        % Add left (negative coefficients)
        for k = 1:numel(leftTerms)
            [coef, met] = parseTerm(leftTerms{k});
            if coef ~= 0
                ii = addMet(met);
                I(end+1) = ii; %#ok<AGROW>
                J(end+1) = j;  %#ok<AGROW>
                V(end+1) = -coef; %#ok<AGROW>
            end
        end
        % Add right (positive coefficients)
        for k = 1:numel(rightTerms)
            [coef, met] = parseTerm(rightTerms{k});
            if coef ~= 0
                ii = addMet(met);
                I(end+1) = ii; %#ok<AGROW>
                J(end+1) = j;  %#ok<AGROW>
                V(end+1) = +coef; %#ok<AGROW>
            end
        end

        rxnIDs(j) = "R" + string(j); % simple IDs R1..Rn
    end

    % Build S
    S = sparse(I, J, V, length(mets), nRxn);

    % -------- 3) Build bounds and metadata  --------
    % COBRA defaults
    LB_DEFAULT_IRREV = 0;       % irreversible lower bound
    LB_DEFAULT_REV   = -1000;   % reversible lower bound
    UB_DEFAULT       = 1000;    % universal upper bound

    lb = zeros(nRxn,1); ub = UB_DEFAULT*ones(nRxn,1);

    % Mark reversibility by arrow type
    for j = 1:nRxn
        if arrows(j) == "<=>"
            lb(j) = LB_DEFAULT_REV;
        else
            lb(j) = LB_DEFAULT_IRREV;
        end
    end

    % ----- TUNE BOUNDS (optional, recommended) -----
    % Option A (safe default): keep wide bounds as above.
    % Option B: narrow bounds around 13C-MFA MinBest ± CI (if present).
    %           This anchors your model to the published fluxes.
    useTightBounds = false;  % change to true if you want tight bounds

    if useTightBounds
        minBest = cell2mat(T.MinBest);
        minLB95 = cell2mat(T.MinLB);
        minUB95 = cell2mat(T.MinUB);
        for j = 1:nRxn
            if ~isnan(minLB95(j)) && ~isnan(minUB95(j))
                lb(j) = minLB95(j);
                ub(j) = minUB95(j);
                % If the reaction was written as irreversible but has lb<0,
                % let it be reversible:
                if lb(j) < 0 && arrows(j) == "->"
                    % keep as is; sign now follows MFA range
                end
            end
        end
    end

    % Objective: try to find a biomass reaction by name
    c = zeros(nRxn,1);
    bioIdx = find(contains(rxnNames, 'Biomass', 'IgnoreCase', true), 1, 'first');
    if ~isempty(bioIdx)
        c(bioIdx) = 1;
    end

    % RHS (mass-balanced internal model, usually zeros)
    b = zeros(length(mets),1);

    % Compose COBRA-style struct
    model = struct();
    model.S = S;
    model.rxns = cellstr(rxnIDs);
    model.rxnNames = cellstr(rxnNames);
    model.mets = cellstr(mets);
    model.lb = lb;
    model.ub = ub;
    model.c = c;
    model.b = b;
    model.description = 'Metabolic model parsed from Au et al. 2014 (Minimal Model flux table)';

    % A handy "rev" vector for COBRA (true if lb < 0)
    model.rev = (model.lb < 0);

    % --- basic sanity checks ---
    fprintf('Built model with %d metabolites and %d reactions.\n', length(mets), nRxn);
    nIrrev = sum(arrows == "->"); nRev = sum(arrows == "<=>");
    fprintf('Irreversible: %d, Reversible: %d\n', nIrrev, nRev);

end

function [coef, met] = localParseTerm(termStr)
    % Parse "[coef] MetName" with optional coefficient.
    % Examples: "Gluc.Ext", "2 Pyr", "0.410 Ala", "  ATP "
    s = strtrim(string(termStr));
    if s == "" || s == "NaN"
        coef = 0; met = "";
        return
    end
    % Remove trailing annotations like "(net)" accidentally left
    s = erase(s, {'(net)','(exch)'});
    % Split by whitespace; first token could be a number
    toks = regexp(s, '^\s*([\-+]?\d*\.?\d+(?:[eE][\-+]?\d+)?)?\s*([A-Za-z0-9\.\_\-]+)\s*$', 'tokens', 'once');
    if isempty(toks)
        % Fall back: if there is a '+' inside, it means tokenization failed earlier
        % Try a simpler split
        parts = strtrim(split(s));
        if numel(parts) == 1
            coef = 1; met = parts(1);
        else
            % If first token is numeric
            f = str2double(parts(1));
            if ~isnan(f)
                coef = f;
                met = parts(2);
            else
                coef = 1;
                met = parts(1);
            end
        end
        return
    end
    numTok = toks{1};
    metTok = toks{2};
    if isempty(numTok)
        coef = 1;
    else
        coef = str2double(numTok);
    end
    met = string(metTok);
end
