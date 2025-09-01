function[output] = SubdomProbRounding_ChopOrAdvnpx(A,RoundingStuff,CalculateErrMtrx)

    RoundingSoftware = RoundingStuff{1};    %%% "advanpix" or "chop"
    RoundingType = RoundingStuff{2};        %%% "Mmtrx" or "StieltjessDiagDom"
    RoundingOptions = RoundingStuff{3};     %%% additional parameters for the rounding software and/or the functions
    
    if strcmp(RoundingSoftware,'advanpix')
        NmbDigits = RoundingOptions{1}; mp.Digits(NmbDigits); ScalingOptions = RoundingOptions{2};
        if strcmp(RoundingType,'Mmtrx')
            output = MtrxFactorsLP_Mmtrx_Advanpix(A,ScalingOptions,NmbDigits,CalculateErrMtrx);
        elseif strcmp(RoundingType,'StieltjessDiagDom')
            output = MtrxFactorsLP_Stieltjess_SandwichScaling_Advanpix(A,ScalingOptions,NmbDigits,CalculateErrMtrx);
        end

    elseif strcmp(RoundingSoftware,'chop')
        ChopOptions = RoundingOptions{1}; ScalingOptions = RoundingOptions{2};
        if strcmp(RoundingType,'Mmtrx')
            output = MtrxFactorsLP_Mmtrx_chop(A,ScalingOptions,ChopOptions,CalculateErrMtrx);
        elseif strcmp(RoundingType,'StieltjessDiagDom')
            output = MtrxFactorsLP_StieltjessDiagDom_chop(A,ScalingOptions,ChopOptions,CalculateErrMtrx);
        end
        
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[MyOutput] = MtrxFactorsLP_Mmtrx_Advanpix(A,ScalingOptions,d_s,ReturnErrMtrx)
    if ~issparse(A), disp('error, input not sparse'); end

    %%%%% rescaling
    RescaleToFitRange = ScalingOptions{1}; FractionOfxmaxForScaling = ScalingOptions{3};
    if RescaleToFitRange
        RescaleType = ScalingOptions{2};
        if strcmp(RescaleType,'diag2side')
            [A_Rescaled,D_rowscale,D_colscale] = scale_diag_2side(A); 
        elseif strcmp(RescaleType,'diag2side_symm')
            tol = 1e-4; [A_Rescaled,D_rowscale,D_colscale,~] = scale_diag_2side_symm(A,tol); %%% we should have "D_rowscale = D_colscale"
        end
        xmax = 10^(d_s); mu_FitToRange = xmax/FractionOfxmaxForScaling; 
        A_RescaledToFitRange = mu_FitToRange * A_Rescaled;
    else
        A_RescaledToFitRange = A; D_rowscale = nan; D_colscale = nan; mu_FitToRange = nan;
    end

    [i,j,v] = find(A_RescaledToFitRange); [v_sorted,perm] = sort(v); perm_inv(perm) = 1:numel(perm);
    
    [~,ind_AbsValMin] = min(abs(v_sorted));
    if v_sorted(ind_AbsValMin) > 0 % if we have a positive value being smallest in magn. -> the entry bfr the first occurence (=index "ind_AbsValMin") is the smallest negative entry in abs. val.
        ind_AbsValMin = ind_AbsValMin - 1;
    else % if we have a negative value being smallest in magn. -> if it is repeated, the frst occurence (=index "ind_AbsValMin") is not what we want, we want the last one (so that the next entry is positive)
        while v_sorted(ind_AbsValMin+1) < 0
            ind_AbsValMin = ind_AbsValMin + 1;
        end
    end
    
    % round positive entries
    pwrs_pos = floor(log10(v_sorted(ind_AbsValMin+1:end))) +1; %pwrs_pos = max(pwrs);
    v_mpPos = mp( ceil(v_sorted(ind_AbsValMin+1:end).*10.^(d_s-pwrs_pos)) )./10.^(d_s-pwrs_pos);
    % round negative entries
    pwrs_neg = floor(log10(-v_sorted(1:ind_AbsValMin))) +1; %pwrs_neg = max(pwrs);
    v_mpNeg = mp( ceil(v_sorted(1:ind_AbsValMin).*10.^(d_s-pwrs_neg)) )./10.^(d_s-pwrs_neg);
    v_mp = [v_mpNeg;v_mpPos]; A_rounded = sparse(i,j,v_mp(perm_inv));

    [L,U,p,q] = lu(A_rounded,'vector'); q_inv(q) = 1:length(q); % we need q_inv because the mtrx "Q" should be applied from the left, i.e., to permute columns -- but then we apply it as a permutation of vector, i.e., we should take inverse

    % %%% check 
    DoCheck = false;
    if ReturnErrMtrx
       if ~DoCheck, ErrMtrx_mp = double(A_rounded - A_RescaledToFitRange); else
       ErrMtrx_mp = double(A_rounded - A_RescaledToFitRange); disp(append('norm of Ei_mp_newway < 0: ',num2str(norm(ErrMtrx_mp(ErrMtrx_mp<=0),'fro')))); end
    else, ErrMtrx_mp=nan; end

    MyOutput = {mu_FitToRange, D_rowscale,D_colscale, L, U, p,q_inv, ErrMtrx_mp,A_rounded,A_RescaledToFitRange};
end


function[MyOutput] = MtrxFactorsLP_Stieltjess_SandwichScaling_Advanpix(A,ScalingOptions,d_s,ReturnErrMtrx)
    if ~issparse(A), disp('error, input not sparse'); end

    %%%%% rescaling & rounding
    RescaleToFitRange = ScalingOptions{1}; RescaleType = ScalingOptions{2}; FractionOfxmaxForScaling = ScalingOptions{3};
    TreatDiagExplicitly = ScalingOptions{4};
    xmax = 10^(d_s); mu_FitToRange = xmax/FractionOfxmaxForScaling;
    if TreatDiagExplicitly
        %%% fit to range - put zeros on the diagonal and at the end 
        D_rowscale = spdiags( 1./sqrt(spdiags(A,0)) , 0 , size(A,1),size(A,2) ); D_colscale = D_rowscale;
        A_Rescaled_wNaturalDiag = D_rowscale * A * D_colscale; A_Rescaled_woDiag = spdiags(zeros(size(A,1),1),0,A_Rescaled_wNaturalDiag);
        if RescaleToFitRange
            A_RescaledToFitRange_woDiag = mu_FitToRange * A_Rescaled_woDiag;
            if ReturnErrMtrx, A_RescaledToFitRange_wDiag = spdiags(mu_FitToRange*ones(size(A,1),1),0,mu_FitToRange*A_Rescaled_wNaturalDiag); else, A_RescaledToFitRange_wDiag = nan; end
        else
            A_RescaledToFitRange_woDiag = A_Rescaled_woDiag;
            if ReturnErrMtrx, A_RescaledToFitRange_wDiag = spdiags(mu_FitToRange*ones(size(A,1),1),0,mu_FitToRange*A_Rescaled_wNaturalDiag); else, A_RescaledToFitRange_wDiag = nan; end
        end

        %%% round
        [i,j,v] = find(A_RescaledToFitRange_woDiag);
        % round off-diag entries (they're all negative)
        pwrs_neg = floor(log10(-v)) +1; % pwrs_neg = max(pwrs);
        v_mp = mp( ceil(v.*10.^(d_s-pwrs_neg)) )./10.^(d_s-pwrs_neg);
        A_Rounded_bfrAddingDiag = sparse(i,j,v_mp); 
        if RescaleToFitRange
            A_Rounded_mp = mp( spdiags(mu_FitToRange*ones(size(A,1),1),0,A_Rounded_bfrAddingDiag) );
        else
            A_Rounded_mp = mp( spdiags(ones(size(A,1),1),0,A_Rounded_bfrAddingDiag) );
        end
    end
    %%% Notice that we chose "mu_FitToRange" as "xmax" divided by "a power of two", hence there is no chance for rounding error.
    %%% Thereby once we round "A_RescaledToFitRange", we don't need to touch the diagonal and only need to rescale the off-diagonal.  
    %%% Easy to check with the code " norm( spdiags(A_RescaledToFitRange,0) - chop_dense(spdiags(A_RescaledToFitRange,0)) ) "

    [L,U,p,q] = lu(A_Rounded_mp,'vector'); q_inv(q) = 1:length(q); % we need q_inv because the mtrx "Q" should be applied from the left, i.e., to permute columns -- but then we apply it as a permutation of vector, i.e., we should take inverse

    % %%% check 
    if ReturnErrMtrx, ErrMtrx_mp = double(A_Rounded_mp - A_RescaledToFitRange_wDiag); else, ErrMtrx_mp=nan; end

    MyOutput = {mu_FitToRange, D_rowscale,D_colscale, L, U, p,q_inv, ErrMtrx_mp,A_Rounded_mp,A_RescaledToFitRange_wDiag};
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[MyOutput] = MtrxFactorsLP_Mmtrx_chop(A,ScalingOptions,ChopOptions,ReturnAddRoundingOutputs)
    if ~issparse(A), disp('error, input not sparse'); end
    
    %%%%% rescaling
    RescaleToFitRange = ScalingOptions{1}; FractionOfxmaxForScaling = ScalingOptions{3};
    if RescaleToFitRange
        RescaleType = ScalingOptions{2};
        if strcmp(RescaleType,'diag2side')
            [A_Rescaled,D_rowscale,D_colscale] = scale_diag_2side(A); 
        elseif strcmp(RescaleType,'diag2side_symm')
            tol = 1e-4; [A_Rescaled,D_rowscale,D_colscale,~] = scale_diag_2side_symm(A,tol); %%% we should have "D_rowscale = D_colscale"
        end
        [~,~,~,xmax,~,~,~] = float_params(ChopOptions.format); 
        mu_FitToRange = xmax/FractionOfxmaxForScaling; A_RescaledToFitRange = mu_FitToRange * A_Rescaled;
    else
        A_RescaledToFitRange = A; D_rowscale = nan; D_colscale = nan; mu_FitToRange = nan;
    end


    %%%%% rounding
    [i,j,v] = find(A_RescaledToFitRange); [v_sorted,perm] = sort(v); perm_inv(perm) = 1:numel(perm);
    
    [~,ind_AbsValMin] = min(abs(v_sorted));
    if v_sorted(ind_AbsValMin) > 0 % if we have a positive value being smallest in magn. -> the entry bfr the first occurence (=index "ind_AbsValMin") is the smallest negative entry in abs. val.
        ind_AbsValMin = ind_AbsValMin - 1;
    else % if we have a negative value being smallest in magn. -> if it is repeated, the frst occurence (=index "ind_AbsValMin") is not what we want, we want the last one (so that the next entry is positive)
        while v_sorted(ind_AbsValMin+1) < 0
            ind_AbsValMin = ind_AbsValMin + 1;
        end
    end

    % round positive entries -> towards +inf
    ChopOptions.round = 2; chop_dense([],ChopOptions)
    v_PosRounded = chop_dense(v_sorted(ind_AbsValMin+1:end));
    % round negative entries -> towards 0
    ChopOptions.round = 4; chop_dense([],ChopOptions)
    v_NegRounded = chop_dense(v_sorted(1:ind_AbsValMin));
    v_Rounded = [v_NegRounded;v_PosRounded]; A_rounded = sparse(i,j,v_Rounded(perm_inv));



    %%%%% LU factorization
    % find the permutations for stability and fill-in reduction
    [~,~,p,q] = lu(A_RescaledToFitRange,'vector'); q_inv(q) = 1:length(q); pivoting = 'none';
    % set the chop_sparse (and return rounding to default)
    ChopOptions.round = 1; chop_sparse([],ChopOptions)
    % calculate the LU
    [L,U] = lu_sparseMmtrx_chop(A_rounded(p,q),pivoting,ChopOptions);
    % save the LU by columns (by testing this was fastest for sparse_chop-solves with the factors)
    L_struct = cell(size(A,2),1); U_struct = cell(size(A,2),1);
    for col = 1:size(A,2)
        L_struct{col} = find( L(col+1:end,col) );
        U_struct{col} = find( U(1:col-1,col) );
    end
    
    %%%%% check
    if ReturnAddRoundingOutputs
        ErrMtrx = A_rounded - A_RescaledToFitRange; 
        %if ~isempty(ErrMtrx(ErrMtrx<0)), disp(append('norm of Ei < 0: ',num2str(norm(ErrMtrx(ErrMtrx<=0),'fro')))); end
    else, ErrMtrx=nan; end

    %%%%% output packaging
    MyOutput = {mu_FitToRange, D_rowscale,D_colscale, L,L_struct, U,U_struct, p,q_inv, ErrMtrx,A_rounded,A_RescaledToFitRange};
end



function[MyOutput] = MtrxFactorsLP_StieltjessDiagDom_chop(A,ScalingOptions,ChopOptions,ReturnAddRoundingOutputs)
    if ~issparse(A), disp('error, input not sparse'); end

    %%%%% rescaling & rounding
    RescaleToFitRange = ScalingOptions{1}; RescaleType = ScalingOptions{2}; FractionOfxmaxForScaling = ScalingOptions{3};
    TreatDiagExplicitly = ScalingOptions{4};
    [~,~,~,xmax,~,~,~] = float_params(ChopOptions.format); mu_FitToRange = xmax/FractionOfxmaxForScaling;
    if TreatDiagExplicitly
        %%% fit to range - put zeros on the diagonal and at the end 
        D_rowscale = spdiags( 1./sqrt(spdiags(A,0)) , 0 , size(A,1),size(A,2) ); D_colscale = D_rowscale;
        A_Rescaled_wNaturalDiag = D_rowscale * A * D_colscale; A_Rescaled_woDiag = spdiags(zeros(size(A,1),1),0,A_Rescaled_wNaturalDiag);
        if RescaleToFitRange
            A_RescaledToFitRange_woDiag = mu_FitToRange * A_Rescaled_woDiag;
            if ReturnAddRoundingOutputs, A_RescaledToFitRange_wDiag = spdiags(mu_FitToRange*ones(size(A,1),1),0,mu_FitToRange*A_Rescaled_wNaturalDiag); else, A_RescaledToFitRange_wDiag = nan; end
        else
            A_RescaledToFitRange_woDiag = A_Rescaled_woDiag;
            if ReturnAddRoundingOutputs, A_RescaledToFitRange_wDiag = spdiags(mu_FitToRange*ones(size(A,1),1),0,mu_FitToRange*A_Rescaled_wNaturalDiag); else, A_RescaledToFitRange_wDiag = nan; end
        end

        %%% round
        [i,j,v] = find(A_RescaledToFitRange_woDiag);
        % round off-diag entries towards 0 (they're all negative)
        ChopOptions.round = 4; chop_dense([],ChopOptions)
        v_Rounded = chop_dense(v); A_Rounded_bfrAddingDiag = sparse(i,j,v_Rounded);
        if RescaleToFitRange
            A_Rounded = spdiags(mu_FitToRange*ones(size(A,1),1),0,A_Rounded_bfrAddingDiag);
        else
            A_Rounded = spdiags(ones(size(A,1),1),0,A_Rounded_bfrAddingDiag);
        end
    end
    %%% Notice that we chose "mu_FitToRange" as "xmax" divided by "a power of two", hence there is no chance for rounding error.
    %%% Thereby once we round "A_RescaledToFitRange", we don't need to touch the diagonal and only need to rescale the off-diagonal.  
    %%% Easy to check with the code " norm( spdiags(A_RescaledToFitRange,0) - chop_dense(spdiags(A_RescaledToFitRange,0)) ) "


    %%%%% LU factorization
    [~,~,p,q] = lu(A_Rounded,'vector'); q_inv(q) = 1:length(q); ChosenPrec = ChopOptions.format; pivoting = 'none'; % we need q_inv because the mtrx "Q" should be applied from the left, i.e., to permute columns -- but then we apply it as a permutation of vector, i.e., we should take inverse
    [L,U] = lu_sparseMmtrx_chop(A_Rounded(p,q),pivoting,ChosenPrec);
    L_struct = cell(size(A,2),1); U_struct = cell(size(A,2),1);
    for col = 1:size(A,2)
        L_struct{col} = find( L(col+1:end,col) );
        U_struct{col} = find( U(1:col-1,col) );
    end

    
    %%%%% check
    if ReturnAddRoundingOutputs
        ErrMtrx = A_Rounded - A_RescaledToFitRange_wDiag; 
        %if ~isempty(ErrMtrx(ErrMtrx<0)), disp(append('norm of Ei < 0: ',num2str(norm(ErrMtrx(ErrMtrx<=0),'fro')))); end
    else, ErrMtrx=nan; end

    %%%%% output packaging
    MyOutput = {mu_FitToRange, D_rowscale,D_colscale, L,L_struct, U,U_struct, p,q_inv, ErrMtrx,A_Rounded,A_RescaledToFitRange_wDiag};

end
