% SOMA PARETO
function [Best , FEs , Mig] = SOMA_PARETO(Info , SOMApara , CostFunction)
    % -------------- Extract Information ----------------------------------
    Max_FEs           = Info.Max_FEs;
    dimension         = Info.dimension;
    % -------------- The search range of the function ---------------------
    VarMin            = Info.Search_Range(1);
    VarMax            = Info.Search_Range(2);
    % -------------- Initial Parameters of SOMA ---------------------------
    PopSize           = SOMApara.PopSize;
    N_jump            = SOMApara.N_jump;
    % --------------------- Create Initial Population ---------------------
    pop               = VarMin + rand(dimension,PopSize)*(VarMax - VarMin);
    fit               = CostFunction(pop)';
    FEs               = PopSize;
    [global_cost, id] = min(fit);
    global_pos        = pop(:,id);
    % ---------------- SOMA MIGRATIONS ------------------------------------
    C           = ceil(PopSize * 0.04);
    A           = ceil(PopSize * 0.20);
    D           = ceil(PopSize * 0.16);
    Mig         = 0;
    while (FEs+N_jump <= Max_FEs)
        Mig     = Mig  + 1;
		% ------------ Control parameters ---------------------------------
        PRT     = 0.50 + 0.45*cos(pi*FEs/Max_FEs+pi);
        Step    = 0.35 + 0.15*cos(pi*FEs/Max_FEs);
        % ------------ Sort POP -------------------------------------------
        pop_sort     = [fit ; pop];
        pop_sort     = sortrows(pop_sort')';
        fit          = pop_sort(1,:);
        pop          = pop_sort(2:end,:);
        % ------------ Migrant selection ----------------------------------
        Migrant_idx  = randi([A+1 A+D]);
        Migrant      = pop(:,Migrant_idx);
        % ------------ Leader selection -----------------------------------
        idx_tar      = randi([1 C]);
        Leader       = pop(:,idx_tar);
        % ------------- Moving process ------------------------------------
        offs_path    = [];
        for move     = 1 : N_jump
            k        = move * Step;
            %----- SOMA Mutation ------------------------------------------
            PRTVector = rand(dimension,1) < PRT;
            %PRTVector = (PRTVector-1)*(1-FEs/Max_FEs) + 1;
            %----- SOMA Crossover -----------------------------------------
            offspring =  Migrant + (Leader - Migrant)*k.*PRTVector;
            offs_path = [offs_path   offspring];
        end % END JUMPING
        %----- Checking Boundary and Replaced Outsize Individuals ---------
        number_offs_path = size(offs_path,2);
        for cl = 1 : number_offs_path
            for rw = 1 : dimension
                if  (offs_path(rw,cl) < VarMin) || (offs_path(rw,cl) > VarMax)
                    offs_path(rw,cl)  = VarMin  +  rand*(VarMax - VarMin);
                end
            end
        end
        %----- SOMA Re-Evaluate Fitness Fuction ---------------------------
        new_cost               = CostFunction(offs_path);
        FEs                    = FEs + number_offs_path;
        %----- SOMA Accepting: Place Best Individual to Population --------
        [Migrant_best, idzz]   = min(new_cost);
        %------------------------------------------------------------------
        if  Migrant_best       <= fit(Migrant_idx)
            fit(Migrant_idx)    = Migrant_best;
            pop(:,Migrant_idx)  = offs_path(:,idzz);
            %----- SOMA Update Global_Leader ------------------------------
            if  Migrant_best   <= global_cost
                global_cost     = Migrant_best;
                global_pos      = offs_path(:,idzz);
            end
        end
    end   % END MIGRATIONS (While Loop)
    Best.Value      = global_cost;
    Best.Positon    = global_pos;
end
%--------------------------------------------------------------------------
% This algorithm is programmed according to the descriptions in the papers 
% listed below:

% Link of paper: https://mendel-journal.org/index.php/mendel/article/view/87
% Diep, Q. B., Zelinka, I., & Das, S. (2019, June). Self-organizing migrating algorithm pareto. In Mendel (Vol. 25, No. 1, pp. 111-120).
% DOI:https://doi.org/10.13164/mendel.2019.1.111

% The control parameters PopSize, N_jump, and Step are closely related 
% and greatly affect the performance of the algorithm. Please refer to the 
% above paper to use the correct control parameters.