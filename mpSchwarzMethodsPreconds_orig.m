% clear; clc; close('all');
% 
% %%% choose set-up
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SM_type = 'MS'; SM_nmbsubdoms_PwrOfTwo = 1; SM_nmbiter = 50; SM_relresacc = 1e-12; 
% RM_nmbdigits = 4;
% Advanpix = false; CalcErrMtrx = false; SandwichScaling = false;
% dampingTheta = 1/(SM_nmbsubdoms_PwrOfTwo+1);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nmb_int_gridcols = 50; % number of "unknown" gridcolumns - excluding the left-most and right-most where we have Dirichlet BCs 
% h=1/(nmb_int_gridcols+1); % mesh size
% 
% % %%% do symmetric AdvecDiff based on SzyldFrommer eqn (15,18) p.648
% % G=numgrid('S',nmb_int_gridcols+2);
% % eta = @(x,y) 0.*x+0.*y; alpha = @(x,y) 1+9.*(x+y); beta = @(x,y) 1+9.*(x+y); nu = @(x,y) 0.*x+0.*y; mu = @(x,y) 0.*x+0.*y;
% % A=ReactAdvDiff_Sqr_FD('SzyldFrommer',G,{eta,alpha,beta,nu,mu});
% % rhs = ones(size(A,1),1); u_init = zeros(size(A,1),1);
% % u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35;
% % if SandwichScaling, RM_type = 'Stieltjess_RoundSandwichScaledBi_Facts'; else, RM_type = 'Stieltjess_RoundBi_Facts'; end
% 
% %%% do non-symmetric AdvecDiff based on SzyldGander Fig 2.1
% G=numgrid('S',nmb_int_gridcols+2);
% eta=inline('x.^2.*cos(x+y).^2','x','y'); a=inline('(x+y).^2.*exp(x-y)','x','y'); b1=inline('(y-0.5)','x','y'); b2=inline('-(x-0.5)','x','y');
% A=ReactAdvDiff_Sqr_FD('SzyldGander',G,{eta,a,b1,b2});
% rhs = ones(size(A,1),1); u_init = zeros(size(A,1),1);
% u_plot = zeros(nmb_int_gridcols+2); angle1 = 224; angle2 = 35; RM_type = 'Mmtrx';
% 
% %%% run test
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if strcmp(SM_type,'dAS'), SchwarzMethod = {SM_type, 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc, dampingTheta}; else, SchwarzMethod = {SM_type, 2^SM_nmbsubdoms_PwrOfTwo, SM_nmbiter, SM_relresacc}; end
% RoundingMethod = {RM_type, RM_nmbdigits, Advanpix, CalcErrMtrx}; debug = 0; v = ones(size(rhs));
% 
% [MyPrecondPackage] = mpSchwarzMethodsPrecondsConstruct(A, SchwarzMethod,RoundingMethod, debug);
% [t,precvec] = mpSchwarzMethodsPreconds_test(v, {A, MyPrecondPackage});






function[t_out,Prec_v] = mpSchwarzMethodsPreconds(v, AdditionalArgs)

A = AdditionalArgs{1}; mpSchwarzPrecond = AdditionalArgs{2};
SchwarzMethod = mpSchwarzPrecond{1}; RoundingMethod = mpSchwarzPrecond{2}; MySubdomFactors_LP = mpSchwarzPrecond{3}; MySubdomIndices = mpSchwarzPrecond{4};
SM_type = SchwarzMethod{1}; SM_nmbsubdom_input = SchwarzMethod{2}; if strcmp(SM_type,'dAS'), damping_theta = SchwarzMethod{end}; end %; SM_nmbiter = SchwarzMethod{3}; SM_relresacc = SchwarzMethod{4}


RM_type = RoundingMethod{1}; RM_nmbdigits = RoundingMethod{2}; RM_Advanpix = RoundingMethod{3}; %RM_CalcErrMtrx = RoundingMethod{4};

SM_nmbsubdom_PwrOfTwo = floor(log2(SM_nmbsubdom_input)); SM_nmbsubdom = 2^SM_nmbsubdom_PwrOfTwo; 
if SM_nmbsubdom ~= SM_nmbsubdom_input, disp(append(append(append('The nmb of subdoms needs to be 2^k for some k. Hence we changed the given',num2str(SM_nmbsubdom_input)),' to '),num2str(SM_nmbsubdom))); end
N = length(v); Prec_v = zeros(N,1);
frst_rows = MySubdomIndices{1}; last_rows = MySubdomIndices{2}; 




% %%% calculate M_{dAS}^{-1} * v. See [SzyldFrommer, AlgConvThryRAS,Sec. 3, 2001] &&& [SaadBook, 14.3.3]
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(SM_type,'dAS')
    rhs = v;
    TimerSM = tic;
    for ind_subdom = 1:SM_nmbsubdom    
        %%% do subdomain solve for restricted rhs
        curr_rhs_bfrperm = rhs(frst_rows(ind_subdom):last_rows(ind_subdom));
        
        if strcmp(RM_type,'Mmtrx')
            curr_Lfact = MySubdomFactors_LP{ind_subdom}{1}; curr_Ufact = MySubdomFactors_LP{ind_subdom}{2}; curr_Pvec = MySubdomFactors_LP{ind_subdom}{3}; curr_Qvec_inv = MySubdomFactors_LP{ind_subdom}{4};
            curr_rhs = curr_rhs_bfrperm(curr_Pvec); 
            if RM_Advanpix
                mp.Digits(RM_nmbdigits); curr_rhs_mp = mp(curr_rhs);
                sol_perm = curr_Ufact \ ( curr_Lfact \ curr_rhs_mp); curr_u  = double(sol_perm(curr_Qvec_inv));
            else
                sol_perm = curr_Ufact \ ( curr_Lfact \ curr_rhs); curr_u  = double(sol_perm(curr_Qvec_inv));
            end
            % %%% check
            % Ai = A(frst_rows(ind_subdom):last_rows(ind_subdom),frst_cols(ind_subdom):last_cols(ind_subdom));
            % curr_Qvec(curr_Qvec_inv) = 1:length(curr_Qvec_inv); %q_inv(q) = 1:length(q);
            % norm(  Ai(curr_Pvec,curr_Qvec) - curr_Lfact*curr_Ufact,'fro') / norm(Ai,'fro')
            % curr_u_check = Ai \ curr_rhs_bfrperm; norm(curr_u-curr_u_check)/norm(curr_u_check) * 1/condest(Ai)
        elseif strcmp(RM_type,'Stieltjess_RoundSandwichScaledBi_Facts')
            curr_Lfact = MySubdomFactors_LP{ind_subdom}{1}; curr_Ufact = MySubdomFactors_LP{ind_subdom}{2}; curr_Pvec = MySubdomFactors_LP{ind_subdom}{3}; curr_Qvec_inv = MySubdomFactors_LP{ind_subdom}{4};
            Ai_diag = MySubdomFactors_LP{ind_subdom}{5}; curr_rhs_bfrperm_scaled = 1 ./ sqrt(Ai_diag) .* curr_rhs_bfrperm;
            curr_rhs = curr_rhs_bfrperm_scaled(curr_Pvec); 
            if RM_Advanpix
                warning('off','all')
                mp.Digits(RM_nmbdigits); curr_rhs_mp = mp(curr_rhs); sol_perm = curr_Ufact \ ( curr_Lfact \ curr_rhs_mp); 
                warning('on','all')
            else
                sol_perm = curr_Ufact \ ( curr_Lfact \ curr_rhs);
            end
             curr_u_unscaled = double(sol_perm(curr_Qvec_inv)); curr_u = 1 ./ sqrt(Ai_diag) .* curr_u_unscaled;
            %%% check
            % Ai = A(frst_rows(ind_subdom):last_rows(ind_subdom),frst_cols(ind_subdom):last_cols(ind_subdom));
            % curr_Qvec(curr_Qvec_inv) = 1:length(curr_Qvec_inv); %q_inv(q) = 1:length(q);
            % Ai_sandwichScaled = diag(1./sqrt(Ai_diag)) * Ai * diag(1./sqrt(Ai_diag));
            % norm(  Ai_sandwichScaled(curr_Pvec,curr_Qvec) - curr_Lfact*curr_Ufact,'fro')
            % curr_u_check = Ai \ curr_rhs_bfrperm; norm(curr_u-curr_u_check)
        elseif strcmp(RM_type,'Stieltjess_RoundBi_Facts')
            ... % need to fill in
        elseif strcmp(RM_type,'Stieltjess_RoundBi_NoFacts') % use GMRES for subdomain solves
            ... % need to fill in
        end
    
        %%% prolong the subdomain solution appropriately and add it to the output vector
        if ind_subdom == 1 %%% update my global solution for subdomain 1
            Prec_v = Prec_v + damping_theta * [curr_u; zeros(N-last_rows(ind_subdom),1)];
        elseif ind_subdom == SM_nmbsubdom %%% update my global solution for subdomain N
            Prec_v = Prec_v + damping_theta * [zeros(frst_rows(ind_subdom)-1,1); curr_u];
        else %%% update my global solution for middle subdomains
            Prec_v = Prec_v + damping_theta * [zeros(frst_rows(ind_subdom)-1,1); curr_u; zeros(N-last_rows(ind_subdom),1)];
        end
    end
    t = toc(TimerSM);




% %%% calculate M_{RAS}^{-1} * v. See [SzyldFrommer, AlgConvThryRAS,Sec. 3, 2001] &&& [SaadBook, Sec. 14.3.3]
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(SM_type,'RAS')
    frst_rows_RAS = MySubdomIndices{3}; last_rows_RAS = MySubdomIndices{4}; rhs = v;
    TimerSM = tic;
    for ind_subdom = 1:SM_nmbsubdom
        %%% do subdomain solve for restricted rhs
        curr_rhs_bfrperm = rhs(frst_rows(ind_subdom):last_rows(ind_subdom));
        if strcmp(RM_type,'Mmtrx')
            curr_Lfact = MySubdomFactors_LP{ind_subdom}{1}; curr_Ufact = MySubdomFactors_LP{ind_subdom}{2}; curr_Pvec = MySubdomFactors_LP{ind_subdom}{3}; curr_Qvec_inv = MySubdomFactors_LP{ind_subdom}{4};
            curr_rhs = curr_rhs_bfrperm(curr_Pvec); sol_perm = curr_Ufact \ ( curr_Lfact \ curr_rhs); curr_u  = sol_perm(curr_Qvec_inv);
        elseif strcmp(RM_type,'Stieltjess_RoundBi_Facts')
            curr_Lfact = MySubdomFactors_LP{ind_subdom}{1}; curr_Ufact = MySubdomFactors_LP{ind_subdom}{2}; curr_Pvec = MySubdomFactors_LP{ind_subdom}{3}; curr_Qvec_inv = MySubdomFactors_LP{ind_subdom}{4};
            Ai_diag = MySubdomFactors_LP{ind_subdom}{5}; curr_rhs_bfrperm = 1 ./ Ai_diag .* curr_rhs_bfrperm;
            curr_rhs = curr_rhs_bfrperm(curr_Pvec); sol_perm = curr_Ufact \ ( curr_Lfact \ curr_rhs); curr_u = sol_perm(curr_Qvec_inv);
        elseif strcmp(RM_type,'Stieltjess_RoundSandwichScaledBi_Facts')
            curr_Lfact = MySubdomFactors_LP{ind_subdom}{1}; curr_Ufact = MySubdomFactors_LP{ind_subdom}{2}; curr_Pvec = MySubdomFactors_LP{ind_subdom}{3}; curr_Qvec_inv = MySubdomFactors_LP{ind_subdom}{4};
            Ai_diag = MySubdomFactors_LP{ind_subdom}{5}; curr_rhs_bfrperm = 1 ./ sqrt(Ai_diag) .* curr_rhs_bfrperm;
            curr_rhs = curr_rhs_bfrperm(curr_Pvec); sol_perm = curr_Ufact \ ( curr_Lfact \ curr_rhs); curr_u_unscaled = sol_perm(curr_Qvec_inv); curr_u = 1 ./ sqrt(Ai_diag) .* curr_u_unscaled;
        elseif strcmp(RM_type,'Stieltjess_RoundBi_NoFacts') % use GMRES for subdomain solves
            ... % need to fill in
        end
        % % check
        % Ai = A(frst_rows(ind_subdom):last_rows(ind_subdom),frst_cols(ind_subdom):last_cols(ind_subdom)); curr_u_check = Ai \ res_prev(frst_rows(ind_subdom):last_rows(ind_subdom));
    
    
        %%% prolong the subdomain solution appropriately and add it to the output vector
        if ind_subdom == 1 %%% update my global solution for subdomain 1
            Prec_v = Prec_v + [curr_u(1:last_rows_RAS(ind_subdom)); zeros(N-last_rows_RAS(ind_subdom),1)];
        elseif ind_subdom == SM_nmbsubdom %%% update my global solution for subdomain N
            frst_row_inside_subdom = frst_rows_RAS(ind_subdom) - frst_rows(ind_subdom) +1;
            Prec_v = Prec_v + [zeros(last_rows_RAS(ind_subdom-1),1); curr_u(frst_row_inside_subdom:end)]; 
        else %%% update my global solution for middle subdomains
            frst_row_inside_subdom = frst_rows_RAS(ind_subdom) - frst_rows(ind_subdom) +1; % nmb_of_rows_in_left_overlap = frst_rows_RAS(ind_subdom) - frst_rows(ind_subdom) -> doesn't include the +1 to start after the overlap
            last_row_inside_subdom = last_rows_RAS(ind_subdom) - frst_rows(ind_subdom) +1; % nmb_of_rows_to_right_overlap = last_rows_RAS(ind_subdom) - frst_rows(ind_subdom) -> doesn't include the +1 to start after the overlap
            Prec_v = Prec_v + [zeros(last_rows_RAS(ind_subdom-1),1); curr_u(frst_row_inside_subdom:last_row_inside_subdom); zeros(N-last_rows_RAS(ind_subdom),1)];
        end
    end
    t = toc(TimerSM);











% %%% calculate M_{MS}^{-1} * v. See [SaadBook, Sec. 14.3.2]
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(SM_type,'MS') 

    TimerSM = tic;
    % solve on the first subdomain, i.e., z = T1*v
    ind_subdom = 1; curr_rhs_bfrperm = v(frst_rows(ind_subdom):last_rows(ind_subdom));
    if strcmp(RM_type,'Mmtrx')
        curr_Lfact = MySubdomFactors_LP{ind_subdom}{1}; curr_Ufact = MySubdomFactors_LP{ind_subdom}{2}; curr_Pvec = MySubdomFactors_LP{ind_subdom}{3}; curr_Qvec_inv = MySubdomFactors_LP{ind_subdom}{4};
        curr_rhs = curr_rhs_bfrperm(curr_Pvec); sol_perm = curr_Ufact \ ( curr_Lfact \ curr_rhs); curr_z_subdom1  = sol_perm(curr_Qvec_inv);
    elseif strcmp(RM_type,'Stieltjess_RoundBi_Facts')
        curr_Lfact = MySubdomFactors_LP{ind_subdom}{1}; curr_Ufact = MySubdomFactors_LP{ind_subdom}{2}; curr_Pvec = MySubdomFactors_LP{ind_subdom}{3}; curr_Qvec_inv = MySubdomFactors_LP{ind_subdom}{4};
        Ai_diag = MySubdomFactors_LP{ind_subdom}{5}; curr_rhs_bfrperm = 1 ./ Ai_diag .* curr_rhs_bfrperm;
        curr_rhs = curr_rhs_bfrperm(curr_Pvec); sol_perm = curr_Ufact \ ( curr_Lfact \ curr_rhs); curr_z_subdom1 = sol_perm(curr_Qvec_inv);
    elseif strcmp(RM_type,'Stieltjess_RoundSandwichScaledBi_Facts')
        curr_Lfact = MySubdomFactors_LP{ind_subdom}{1}; curr_Ufact = MySubdomFactors_LP{ind_subdom}{2}; curr_Pvec = MySubdomFactors_LP{ind_subdom}{3}; curr_Qvec_inv = MySubdomFactors_LP{ind_subdom}{4};
        Ai_diag = MySubdomFactors_LP{ind_subdom}{5}; curr_rhs_bfrperm = 1 ./ sqrt(Ai_diag) .* curr_rhs_bfrperm;
        curr_rhs = curr_rhs_bfrperm(curr_Pvec); sol_perm = curr_Ufact \ ( curr_Lfact \ curr_rhs); curr_z_unscaled = sol_perm(curr_Qvec_inv); curr_z_subdom1 = 1 ./ sqrt(Ai_diag) .* curr_z_unscaled;
    elseif strcmp(RM_type,'Stieltjess_RoundBi_NoFacts') % use GMRES for subdomain solves
        ... % need to fill in
    end
    curr_z_prev = [curr_z_subdom1; zeros(N-last_rows(ind_subdom),1)];
    
    
    
    % continue for the remaining subdomains with z = z + Ti*(v-A*z)
    for ind_subdom = 2:SM_nmbsubdom    
    
        %%% do subdomain solve for restricted rhs
        curr_rhs_bfrperm_global = v - A*curr_z_prev; 
        curr_rhs_bfrperm = curr_rhs_bfrperm_global(frst_rows(ind_subdom):last_rows(ind_subdom));
        if strcmp(RM_type,'Mmtrx')
            curr_Lfact = MySubdomFactors_LP{ind_subdom}{1}; curr_Ufact = MySubdomFactors_LP{ind_subdom}{2}; curr_Pvec = MySubdomFactors_LP{ind_subdom}{3}; curr_Qvec_inv = MySubdomFactors_LP{ind_subdom}{4};
            curr_rhs = curr_rhs_bfrperm(curr_Pvec); sol_perm = curr_Ufact \ ( curr_Lfact \ curr_rhs); curr_z_subdom_i  = sol_perm(curr_Qvec_inv);
        elseif strcmp(RM_type,'Stieltjess_RoundBi_Facts')
            curr_Lfact = MySubdomFactors_LP{ind_subdom}{1}; curr_Ufact = MySubdomFactors_LP{ind_subdom}{2}; curr_Pvec = MySubdomFactors_LP{ind_subdom}{3}; curr_Qvec_inv = MySubdomFactors_LP{ind_subdom}{4};
            Ai_diag = MySubdomFactors_LP{ind_subdom}{5}; curr_rhs_bfrperm = 1 ./ Ai_diag .* curr_rhs_bfrperm;
            curr_rhs = curr_rhs_bfrperm(curr_Pvec); sol_perm = curr_Ufact \ ( curr_Lfact \ curr_rhs); curr_z_subdom_i = sol_perm(curr_Qvec_inv);
        elseif strcmp(RM_type,'Stieltjess_RoundSandwichScaledBi_Facts')
            curr_Lfact = MySubdomFactors_LP{ind_subdom}{1}; curr_Ufact = MySubdomFactors_LP{ind_subdom}{2}; curr_Pvec = MySubdomFactors_LP{ind_subdom}{3}; curr_Qvec_inv = MySubdomFactors_LP{ind_subdom}{4};
            Ai_diag = MySubdomFactors_LP{ind_subdom}{5}; curr_rhs_bfrperm = 1 ./ sqrt(Ai_diag) .* curr_rhs_bfrperm;
            curr_rhs = curr_rhs_bfrperm(curr_Pvec); sol_perm = curr_Ufact \ ( curr_Lfact \ curr_rhs); curr_u_unscaled = sol_perm(curr_Qvec_inv); curr_z_subdom_i = 1 ./ sqrt(Ai_diag) .* curr_u_unscaled;
        elseif strcmp(RM_type,'Stieltjess_RoundBi_NoFacts') % use GMRES for subdomain solves
            ... % need to fill in
        end
        % % check
        % Ai = A(frst_rows(ind_subdom):last_rows(ind_subdom),frst_cols(ind_subdom):last_cols(ind_subdom)); curr_u_check = Ai \ res_prev(frst_rows(ind_subdom):last_rows(ind_subdom));
    
        %%% return to the global coordinates
        if ind_subdom == SM_nmbsubdom %%% update my global solution for subdomain N
            curr_z_updt = [zeros(frst_rows(ind_subdom)-1,1); curr_z_subdom_i];
        else %%% update my global solution for middle subdomains
            curr_z_updt = [zeros(frst_rows(ind_subdom)-1,1); curr_z_subdom_i; zeros(N-last_rows(ind_subdom),1)];
        end
    
        %%% update curr_z
        curr_z = curr_z_prev + curr_z_updt; 
        curr_z_prev = curr_z;
    end
    Prec_v = curr_z;
    t = toc(TimerSM);
    
end

t_out = [t;0];
end

