function [AddOutputs,x, res_norms] = MyGMRES_PrecApply(A, b, PrecInfo, Prec, x0, MaxIt, RelResTol )
	% Implementation of the GMRES method by Schultz and Saad, python adaptaion of "Gander_Iter_Methods_lecture_notes" - page 109 (Chapter 3)    
    % 
    % Parameters
    % ----------
    % A : matrix
    %     N-by-N system matrix (nonsingular)
    % b : vector
    %     N-by-1 vector of the right-hand side
    %
    % PrecInfo: structurte given by PrecInfo = {TypeOfPrec,PrecForm,PrecExtraArgs} with the following options :
    %     - TypeOfPrec = 'noprec' or 'R' or 'L' -> use no preconditioning or left or right preconditioniing
    %     - PrecForm = 'Matrix' or 'Handle'   ~ means "M" is either a concrete matrix or a handle of a function "v \mapsto Prec*v"
    %     - PrecExtraArgs = additional arguments for the call of the routine "Prec(v,PrecExtraArgs)"
    %
    % Prec : numpy matrix
    %     the preconditioner or the function "v \mapsto Prec*v"
    % 
    % x0 : N-by-1 vector of the initial guess or 'false' if no initial guess given
    %
    % tol : float
    %     absolute tolerance for residual convergence
    % nmb_iter : int
    %     number of iterations before GMRES is stopped
    % print_out : bool
    %     if True, the message is printed out as this routine finishes        
    % 
    % Returns
    % -------
    % x : numpy vector
    %     the approximate solution to "Ax=b"
    % res_norms : numpy vector
    %     the vector with norms of residuals throughout the "nmb_iter" iterations of GMRES
    % 
    % 
    % Notes
    % -------
    % - A can be both sparse or dense format
    % 
    % - we use Givens rotations for the QR factorization of the upper Hessenberg matrix H_{k+1,k}. As it si upper Hessenberg,
    % each column needs only one rotation, i.e., one cosine and one sine. We store these in vectors "cos, sin"
    % 
    % - we store the Givens rotations (in "cos, sin") because we in each iteration only update the matrix H_{k+1,k} and its
    % R-factor (from QR factoriz.). We don't recompute these from scratch every iteration
    % 
    % - we don't solve the LS problem 
    %     H_{k+1,k} y_k ~ rhs
    % because we can compute the residual of this problem explicitly based on "rhs" - see the lecture notes, page 110 (Chapter 3). 
    % 
    % - as a result, we only need to update the quantities "rhs, Q, H, R" and once the residual is small enough, we onyl then compute
    % the actual approximation "x"


if ischar(PrecInfo{1})
    
    if strcmp(PrecInfo{1},'noprec')
        TimerGMRES = tic;
        [NmbIters, x, res_norms] = nonprec_gmres(A, b, x0, MaxIt, RelResTol, false);
        GMRESTime = toc(TimerGMRES); AddOutputs = {'GMRES total time',GMRESTime,'nmb iter',NmbIters, };
        %[x, res_norms] = nonprec_gmres_Martin(A,b,RelResTol,MaxIt);
    elseif strcmp(PrecInfo{1},'L')
        TimerGMRES = tic;
        [t,x, res_norms] = leftprec_gmres(A, b, Prec, PrecInfo(2:end), x0, MaxIt, RelResTol, false);
        GMRESTime = toc(TimerGMRES); AddOutputs = {'P1 Time',t(1),'P2 Time',t(2),'GMRES total time',GMRESTime,'nmb iter',t(3)};
    elseif strcmp(PrecInfo{1},'R')
        TimerGMRES = tic;
        [t,x, res_norms] = rightprec_gmres(A, b, Prec, PrecInfo(2:end), x0, MaxIt, RelResTol, false);
        GMRESTime = toc(TimerGMRES); AddOutputs = {'P1 Time',t(1),'P2 Time',t(2),'GMRES total time',GMRESTime,'nmb iter',t(3)};
    end
    
end
    
end






%---------------------------------------------------------------------
% Non-preconditioned GMRES
% --------------------------------------------------------------------
function [NmbOfIters, x, res_norms] = nonprec_gmres(A, b, x_0, nmb_iter, tol, print_out)
    % Implementation of the GMRES method by Schultz and Saad, python adaptaion of "Gander_Iter_Methods_lecture_notes" - page 109 (Chapter 3)    
    % 
    % Parameters
    % ----------
    % A : matrix
    %     the system matrix (nonsingular). Can be either sparse format or normal numpy matrix
    % b : numpy vector
    %     the right-hand side
    % tol : float
    %     absolute tolerance for residual convergence
    % nmb_iter : int
    %     number of iterations before GMRES is stopped
    % print_out : bool
    %     if True, the message is printed out as this routine finishes        
    % 
    % Returns
    % -------
    % x : numpy vector
    %     the approximate solution to "Ax=b"
    % res_norms : numpy vector
    %     the vector with norms of residuals throughout the "nmb_iter" iterations of GMRES
    % 
    % 
    % Notes
    % -------
    % - A can be both sparse or dense format
    % 
    % - we use Givens rotations for the QR factorization of the upper Hessenberg matrix H_{k+1,k}. As it si upper Hessenberg,
    % each column needs only one rotation, i.e., one cosine and one sine. We store these in vectors "cos, sin"
    % 
    % - we store the Givens rotations (in "cos, sin") because we in each iteration only update the matrix H_{k+1,k} and its
    % R-factor (from QR factoriz.). We don't recompute these from scratch every iteration
    % 
    % - we don't solve the LS problem 
    %     H_{k+1,k} y_k ~ rhs
    % because we can compute the residual of this problem explicitly based on "rhs" - see the lecture notes, page 110 (Chapter 3). 
    % 
    % - as a result, we only need to update the quantities "rhs, Q, H, R" and once the residual is small enough, we onyl then compute
    % the actual approximation "x"

    
    
    N = length(b); k = 1; iterate_bool = 1; %Initialize the iteration index
   
    if islogical(x_0) % Initialize the vercor of residuals
        res = b;
    else
        res = b - A*x_0; 
    end
    
    res_norms = zeros(nmb_iter,1); curr_res_norm = norm(res,'fro');
    res_norms(1) = curr_res_norm; rhs = res; % Initialize the right-hand side
    
    if curr_res_norm == 0
        x = zeros(N,1); res_norms = zeros(N,1);
        return
    end
    
    RelTol = tol*curr_res_norm;
    Q_previter = res / curr_res_norm; % First column of Q
    rhs_previter = curr_res_norm; % First "rhs" for the Hessenberg matrix problem
    H_previter = eye(1)*42; R_previter = eye(1)*42; % no need to initialize
    

    % we initialize the Givens rotation "sine" and "cosine" vectors. As for every column of H_{k+1,k} there is only one element
    % below the diagonal we have only one rotation per column. 
    % This rotation for "i"-th column in H_{k+1,k} is given by "cos[i], sin[i]"
    cos = zeros(nmb_iter,1); sin = zeros(nmb_iter,1);



    
    while iterate_bool  % GMRES iterations

        % ##########################################################
        % # PART 1: Arnoldi process " A Q_k = Q_{k+1} H_{k+1,k} "
        % ##########################################################
        Q = zeros(N,k+1); Q(:,1:k) = Q_previter; H = zeros(k+1,k); H(1:k,1:k-1) = H_previter;

        q_kp1 = A * Q(:,k);
        
        for j = 1:k   % Gram-Schmidt on " q_kp1 = A*Q(:,k) "
            H(j,k) = dot(q_kp1 , Q(:,j));
            q_kp1 = q_kp1 - H(j,k) * Q(:,j);
        end

        H(k+1,k) = norm( q_kp1, 'fro' ); % Lower-diagonal element
        Q(:,k+1) = q_kp1 / H(k+1,k);  % Compute the vector Q(:,k+1)
        % disp(k)
        % disp(Q)



        % ########################################################################################################
        % # PART 2: Assemble the R-factor of the Hessenberg matrix "H_{k+1,k}" by Givens rotations
        % ########################################################################################################
        R = zeros(k+1,k); R(1:k,1:k-1) = R_previter;
        R(:,k) = H(:,k); % New column of R
        
        
        for j = 1:k-1 % Previous Givens rotations on the last "newly added" column
            Rjk = cos(j)*R(j,k) + sin(j)*R(j+1,k);
            R(j+1,k) = -sin(j)*R(j,k) + cos(j)*R(j+1,k);
            R(j,k) = Rjk;
        end

        R_lastcol = R(:,k); % The last column of R
        
        % Next compute and apply the "new" Givens rotation to zero-out the entry "h_{k+1}" of the last col of "H_{k+1,k}"
        cos(k) = R_lastcol(k)   / sqrt( R_lastcol(k)^2 + R_lastcol(k+1)^2 );
        sin(k) = R_lastcol(k+1) / sqrt( R_lastcol(k)^2 + R_lastcol(k+1)^2 );
        
        % Checking that the Givens rotation are correct
        % R_kp1_k_aux = -sin(k)*R_lastcol(k) + cos(k)*R_lastcol(k+1)
        % disp('Checking the Givens rotation - should be zero : ',R_kp1_k_aux)

        R(k,k) = cos(k)*R_lastcol(k) + sin(k)*R_lastcol(k+1);
        R(k+1,k) = 0;



        
        
        % ##########################################################
        % # PART 3: Update rhs and residuals
        % ##########################################################
        rhs = zeros(k+1,1); rhs(1:k) = rhs_previter;
        rhs(k+1) = -sin(k)*rhs_previter(k);
        rhs(k)   =  cos(k)*rhs_previter(k);
        % disp(rhs_previter); disp(rhs); disp(sin,cos)
        
        res_norms(k+1) = abs( rhs(k+1) ); % Compute the new norm of the residual
        curr_res_norm = res_norms(k+1);




        k = k+1; % Update iteration index
        rhs_previter = rhs; R_previter = R; H_previter = H; Q_previter = Q;

        if curr_res_norm <= RelTol 
            iterate_bool = false;
            if print_out
                disp('GMRES terminated at iteration ',k,' with the residual norm : ',curr_res_norm); disp(40*'-')
            end
        % else:
        %     print('Iteration ',k,' with residual norm ',curr_res_norm,' with tol = ',tol)
        end

        if k > nmb_iter
            iterate_bool = false;
            if print_out
                disp('GMRES terminated at iteration ',k,' with the residual norm : ',curr_res_norm); disp(40*'-')
            end
        end
     
    end
    
   
    y = linsolve(R,rhs); % Compute the coefficients y_j

    if islogical(x_0) % Add the initial guess
        x = Q_previter(:,1:k-1) * y; % Compute x as sum_j y_j q_j - the last "just added" column of Q shouldnt be used
    else
        x = x_0 + Q_previter(:,1:k-1) * y; % Compute x as sum_j y_j q_j - the last "just added" column of Q shouldnt be used            
    end

    NmbOfIters = k-1;
end












































%---------------------------------------------------------------------
% Right-preconditioned GMRES
% --------------------------------------------------------------------
function [timings,x, res_norms] = rightprec_gmres(A, b, M, PrecAddInfo, x_0, nmb_iter, tol, print_out)
    % Implementation of the GMRES method by Schultz and Saad, python adaptaion of "Gander_Iter_Methods_lecture_notes" - page 109 (Chapter 3)
    % But here we add the right preconditioning - see Saad_Iter_Methods_book - page 284 (Section 9.3)
    % 
    % Parameters
    % ----------
    % A : numpy matrix
    %     the system matrix (nonsingular)
    % b : numpy vector
    %     the right-hand side
    % M : numpy matrix
    %     the preconditioner or the function "v \mapsto Prec*v"
    % 
    % PrecAddInfo : {PrecForm,PrecExtraArgs},  following options :
    %     - PrecForm='Matrix'   ~ means "M" is a sparse matrix and there is going to be a direct solve.
    %     - PrecForm='Handle'   ~ means "M" is handle of a function "v \mapsto Prec*v" with the caller M(v,PrecExtraArgs)
    % 
    % tol : float
    %     absolute tolerance for residual convergence
    % nmb_iter : int
    %     number of iterations before GMRES is stopped
    % print_out : bool
    %     if True, the message is printed out as this routine finishes
    % 
    % Returns
    % -------
    % x : numpy vector
    %     the approximate solution to "Ax=b"
    % res : numpy vector
    %     the vector with norms of residuals throughout the "nmb_iter" iterations of GMRES
    % 
    % 
    % Notes
    % -------
    % - we use Givens rotations for the QR factorization of the upper Hessenberg matrix H_{k+1,k}. As it si upper Hessenberg,
    % each column needs only one rotation, i.e., one cosine and one sine. We store these in vectors "cos, sin"
    % 
    % - we store the Givens rotations (in "cos, sin") because we in each iteration only update the matrix H_{k+1,k} and its
    % R-factor (from QR factoriz.). We don't recompute these from scratch every iteration
    % 
    % - we don't solve the LS problem 
    %     H_{k+1,k} y_k ~ rhs
    % because we can compute the residual of this problem explicitly based on "rhs" - see the lecture notes, page 110 (Chapter 3). 
    % 
    % - as a result, we only need to update the quantities "rhs, Q, H, R" and once the residual is small enough, we onyl then compute
    % the actual approximation "x"
    % 
    % - compared to the non-preconditoned version, the left preconditioned version still "operates" with the desired residuals (although implicitly)
    % in contrast to the left preconditioned version. Thus we start with the "normal" residual 
    % 
    % - compared to the non-preconditoned version, we need to change the "system matrix" in the Arnoldi process
    % 
    % - compared to the non-preconditoned version, after the iteration terminated we need to retrieve the solution "x_m" with an additional solve
    % since we will have
    %     x_m = x_0 + M^{-1} * V_m y_m
    % where y_m was obtained from the oslution of the projectred problem with the upper Heessenberg matrix H_{m+1,m}   

    
    
    N = length(b); k = 1; iterate_bool = 1; %Initialize the iteration index
    if strcmp(PrecAddInfo{1},'Matrix')
        PrecForm = 'Matrix';  TimeP1 = 0; TimeP2 = 0;
    elseif strcmp(PrecAddInfo{1},'Handle')
        PrecForm = 'Handle'; PrecExtraArgs = PrecAddInfo{2}; TimeP1 = 0; TimeP2 = 0;
    end

   
    if islogical(x_0) % Initialize the vercor of residuals
        res = b;
    else
        res = b - A*x_0; % the sparse mat-vec is also ".dot()"
    end
    
    res_norms = zeros(nmb_iter,1); curr_res_norm = norm(res,'fro');
    res_norms(1) = curr_res_norm; rhs = res; % Initialize the right-hand side
    
    if curr_res_norm == 0
        x = zeros(N,1); res_norms = zeros(N,1);
        return
    end
        
    RelTol = tol*curr_res_norm;
    Q_previter = res / curr_res_norm; % First column of Q
    rhs_previter = curr_res_norm; % First "rhs" for the Hessenberg matrix problem
    H_previter = eye(1)*42; R_previter = eye(1)*42; % no need to initialize
    

    % we initialize the Givens rotation "sine" and "cosine" vectors. As for every column of H_{k+1,k} there is only one element
    % below the diagonal we have only one rotation per column. 
    % This rotation for "i"-th column in H_{k+1,k} is given by "cos[i], sin[i]"
    cos = zeros(nmb_iter,1); sin = zeros(nmb_iter,1);
    
    
    while iterate_bool  % GMRES iterations
        
        % ##########################################################
        % # PART 1: Arnoldi process " A Q_k = Q_{k+1} H_{k+1,k} "
        % ##########################################################
        Q = zeros(N,k+1); Q(:,1:k) = Q_previter; H = zeros(k+1,k); H(1:k,1:k-1) = H_previter;

        % "q_kp1_aux = np.linalg.solve(M, Q[:,k -1])"
        if strcmp(PrecForm,'Matrix')
            TimerP1 = tic; q_kp1_aux = M*Q(:,k); t = toc(TimerP1);
            TimeP1 = TimeP1 + t;
        elseif strcmp(PrecForm,'Handle')
            [t,q_kp1_aux] = M(Q(:,k),PrecExtraArgs); TimeP1 = TimeP1 + t(1); TimeP2 = TimeP2 + t(2);
        end
        q_kp1 = A * q_kp1_aux;
        

        for j = 1:k   % Gram-Schmidt on " q_kp1 = A*Q(:,k) "
            H(j,k) = dot(q_kp1 , Q(:,j));
            q_kp1 = q_kp1 - H(j,k) * Q(:,j);
        end

        H(k+1,k) = norm( q_kp1, 'fro' ); % Lower-diagonal element
        Q(:,k+1) = q_kp1 / H(k+1,k);  % Compute the vector Q(:,k+1)




        % ########################################################################################################
        % # PART 2: Assemble the R-factor of the Hessenberg matrix "H_{k+1,k}" by Givens rotations
        % ########################################################################################################
        R = zeros(k+1,k); R(1:k,1:k-1) = R_previter;
        R(:,k) = H(:,k); % New column of R
        
        
        for j = 1:k-1 % Previous Givens rotations on the last "newly added" column
            Rjk = cos(j)*R(j,k) + sin(j)*R(j+1,k);
            R(j+1,k) = -sin(j)*R(j,k) + cos(j)*R(j+1,k);
            R(j,k) = Rjk;
        end

        R_lastcol = R(:,k); % The last column of R
        
        % Next compute and apply the "new" Givens rotation to zero-out the entry "h_{k+1}" of the last col of "H_{k+1,k}"
        cos(k) = R_lastcol(k)   / sqrt( R_lastcol(k)^2 + R_lastcol(k+1)^2 );
        sin(k) = R_lastcol(k+1) / sqrt( R_lastcol(k)^2 + R_lastcol(k+1)^2 );
        
        % Checking that the Givens rotation are correct
        % R_kp1_k_aux = -sin(k)*R_lastcol(k) + cos(k)*R_lastcol(k+1)
        % disp('Checking the Givens rotation - should be zero : ',R_kp1_k_aux)

        R(k,k) = cos(k)*R_lastcol(k) + sin(k)*R_lastcol(k+1);
        R(k+1,k) = 0;


        
        
        % ##########################################################
        % # PART 3: Update rhs and residuals
        % ##########################################################
        rhs = zeros(k+1,1); rhs(1:k) = rhs_previter;
        rhs(k+1) = -sin(k)*rhs_previter(k);
        rhs(k)   =  cos(k)*rhs_previter(k);
        % disp(rhs_previter); disp(rhs); disp(sin,cos)
        
        res_norms(k+1) = abs( rhs(k+1) ); % Compute the new norm of the residual
        curr_res_norm = res_norms(k+1);

        k = k+1; % Update iteration index
        rhs_previter = rhs; R_previter = R; H_previter = H; Q_previter = Q;

        if curr_res_norm <= RelTol 
            iterate_bool = false;
            if print_out
                disp('GMRES terminated at iteration ',k,' with the residual norm : ',curr_res_norm); disp(40*'-')
            end
        % else:
        %     print('Iteration ',k,' with residual norm ',curr_res_norm,' with tol = ',tol)
        end

        if k > nmb_iter
            iterate_bool = false;
            if print_out
                disp('GMRES terminated at iteration ',k,' with the residual norm : ',curr_res_norm); disp(40*'-')
            end
        end
     
    end
    
   
    y = linsolve(R,rhs); % Compute the coefficients y_j
    u = Q_previter(:,1:k-1) * y; % Compute u as sum_j y_j q_j - the last "just added" column of Q shouldnt be used
    if strcmp(PrecForm,'Matrix')
        TimerP1 = tic; x = M*u; t = toc(TimerP1);
        TimeP1 = TimeP1 + t;
    elseif strcmp(PrecForm,'Handle')
        [t,x] = M(u,PrecExtraArgs); TimeP1 = TimeP1 + t(1); TimeP2 = TimeP2 + t(2);
    end

    if ~islogical(x_0) % go back from the preconditioned space and add the initial guess
        x = x0 + x;
    end

timings = [TimeP1,TimeP2,k-1];
end



















%---------------------------------------------------------------------
% Left-preconditioned GMRES
% --------------------------------------------------------------------
function [timings,x, res_norms] = leftprec_gmres(A, b, M, PrecAddInfo, x_0, nmb_iter, tol, print_out)
    % Implementation of the GMRES method by Schultz and Saad, python adaptaion of "Gander_Iter_Methods_lecture_notes" - page 109 (Chapter 3)    
    % 
    % Parameters
    % ----------
    % A : matrix
    %     the system matrix (nonsingular). Can be either sparse format or normal numpy matrix
    % b : numpy vector
    %     the right-hand side
    % M : numpy matrix
    %     the preconditioner or the function "v \mapsto Prec*v"
    % 
    % PrecAddInfo : {PrecForm,PrecExtraArgs},  following options :
    %     - PrecForm='Matrix'   ~ means "M" is a sparse matrix and there is going to be a direct solve.
    %     - PrecForm='Handle'   ~ means "M" is handle of a function "v \mapsto Prec*v" with the caller M(v,PrecExtraArgs)
    % 
    % tol : float
    %     absolute tolerance for residual convergence
    % nmb_iter : int
    %     number of iterations before GMRES is stopped
    % print_out : bool
    %     if True, the message is printed out as this routine finishes        
    % 
    % Returns
    % -------
    % x : numpy vector
    %     the approximate solution to "Ax=b"
    % res : numpy vector
    %     the vector with norms of residuals throughout the "nmb_iter" iterations of GMRES
    % 
    % 
    % Notes
    % -------
    % - A can be both sparse or dense format
    % 
    % - we use Givens rotations for the QR factorization of the upper Hessenberg matrix H_{k+1,k}. As it si upper Hessenberg,
    % each column needs only one rotation, i.e., one cosine and one sine. We store these in vectors "cos, sin"
    % 
    % - we store the Givens rotations (in "cos, sin") because we in each iteration only update the matrix H_{k+1,k} and its
    % R-factor (from QR factoriz.). We don't recompute these from scratch every iteration
    % 
    % - we don't solve the LS problem 
    %     H_{k+1,k} y_k ~ rhs
    % because we can compute the residual of this problem explicitly based on "rhs" - see the lecture notes, page 110 (Chapter 3). 
    % 
    % - as a result, we only need to update the quantities "rhs, Q, H, R" and once the residual is small enough, we onyl then compute
    % the actual approximation "x"

    
    
    N = length(b); k = 1; iterate_bool = 1; %Initialize the iteration index
    if strcmp(PrecAddInfo{1},'Matrix')
        PrecForm = 'Matrix'; TimeP1 = 0; TimeP2 = 0;
    elseif strcmp(PrecAddInfo{1},'Handle')
        PrecForm = 'Handle'; PrecExtraArgs = PrecAddInfo{2}; TimeP1 = 0; TimeP2 = 0;
    end
   
    if islogical(x_0) % Initialize the vercor of residuals
        res_aux = b;
    else
        res_aux = b - A*x_0; % the sparse mat-vec is also ".dot()"
    end

     if strcmp(PrecForm,'Matrix')
        TimerP1 = tic; res = M*res_aux; t = toc(TimerP1);
        TimeP1 = TimeP1 + t;
     elseif strcmp(PrecForm,'Handle')
        [t,res] = M(res_aux,PrecExtraArgs); TimeP1 = TimeP1 + t(1); TimeP2 = TimeP2 + t(2);
     end
    
    res_norms = zeros(nmb_iter,1); curr_res_norm = norm(res,'fro');
    res_norms(1) = curr_res_norm; rhs = res; % Initialize the right-hand side
    
    if curr_res_norm == 0
        x = zeros(N,1); res_norms = zeros(N,1);
        return
    end
        
    RelTol = tol*curr_res_norm;
    Q_previter = res / curr_res_norm; % First column of Q
    rhs_previter = curr_res_norm; % First "rhs" for the Hessenberg matrix problem
    H_previter = eye(1)*42; R_previter = eye(1)*42; % no need to initialize
    

    % we initialize the Givens rotation "sine" and "cosine" vectors. As for every column of H_{k+1,k} there is only one element
    % below the diagonal we have only one rotation per column. 
    % This rotation for "i"-th column in H_{k+1,k} is given by "cos[i], sin[i]"
    cos = zeros(nmb_iter,1); sin = zeros(nmb_iter,1);


    
    while iterate_bool  % GMRES iterations
        
        % ##########################################################
        % # PART 1: Arnoldi process " A Q_k = Q_{k+1} H_{k+1,k} "
        % ##########################################################
        Q = zeros(N,k+1); Q(:,1:k) = Q_previter; H = zeros(k+1,k); H(1:k,1:k-1) = H_previter;

        q_kp1_aux = A * Q(:,k);
        
        % "q_kp1 = np.linalg.solve(M, q_kp1_aux)"
        if strcmp(PrecForm,'Matrix')
            TimerP1 = tic; q_kp1 = M*q_kp1_aux; t = toc(TimerP1);
            TimeP1 = TimeP1 + t;
        elseif strcmp(PrecForm,'Handle')
            [t,q_kp1] = M(q_kp1_aux,PrecExtraArgs); TimeP1 = TimeP1 + t(1); TimeP2 = TimeP2 + t(2);
        end
        
        for j = 1:k   % Gram-Schmidt on " q_kp1 = A*Q(:,k) "
            H(j,k) = dot(q_kp1 , Q(:,j));
            q_kp1 = q_kp1 - H(j,k) * Q(:,j);
        end

        H(k+1,k) = norm( q_kp1, 'fro' ); % Lower-diagonal element
        Q(:,k+1) = q_kp1 / H(k+1,k);  % Compute the vector Q(:,k+1)






        % ########################################################################################################
        % # PART 2: Assemble the R-factor of the Hessenberg matrix "H_{k+1,k}" by Givens rotations
        % ########################################################################################################
        R = zeros(k+1,k); R(1:k,1:k-1) = R_previter;
        R(:,k) = H(:,k); % New column of R
        
        
        for j = 1:k-1 % Previous Givens rotations on the last "newly added" column
            Rjk = cos(j)*R(j,k) + sin(j)*R(j+1,k);
            R(j+1,k) = -sin(j)*R(j,k) + cos(j)*R(j+1,k);
            R(j,k) = Rjk;
        end

        R_lastcol = R(:,k); % The last column of R
        
        % Next compute and apply the "new" Givens rotation to zero-out the entry "h_{k+1}" of the last col of "H_{k+1,k}"
        cos(k) = R_lastcol(k)   / sqrt( R_lastcol(k)^2 + R_lastcol(k+1)^2 );
        sin(k) = R_lastcol(k+1) / sqrt( R_lastcol(k)^2 + R_lastcol(k+1)^2 );
        
        % Checking that the Givens rotation are correct
        % R_kp1_k_aux = -sin(k)*R_lastcol(k) + cos(k)*R_lastcol(k+1)
        % disp('Checking the Givens rotation - should be zero : ',R_kp1_k_aux)

        R(k,k) = cos(k)*R_lastcol(k) + sin(k)*R_lastcol(k+1);
        R(k+1,k) = 0;


        
        
        % ##########################################################
        % # PART 3: Update rhs and residuals
        % ##########################################################
        rhs = zeros(k+1,1); rhs(1:k) = rhs_previter;
        rhs(k+1) = -sin(k)*rhs_previter(k);
        rhs(k)   =  cos(k)*rhs_previter(k);
        % disp(rhs_previter); disp(rhs); disp(sin,cos)
        
        res_norms(k+1) = abs( rhs(k+1) ); % Compute the new norm of the residual
        curr_res_norm = res_norms(k+1);



        k = k+1; % Update iteration index
        rhs_previter = rhs; R_previter = R; H_previter = H; Q_previter = Q;

        if curr_res_norm <= RelTol 
            iterate_bool = false;
            if print_out
                disp('GMRES terminated at iteration ',k,' with the residual norm : ',curr_res_norm); disp(40*'-')
            end
        % else:
        %     print('Iteration ',k,' with residual norm ',curr_res_norm,' with tol = ',tol)
        end

        if k > nmb_iter
            iterate_bool = false;
            if print_out
                disp('GMRES terminated at iteration ',k,' with the residual norm : ',curr_res_norm); disp(40*'-')
            end
        end
     
    end
    
   
    y = linsolve(R,rhs); % Compute the coefficients y_j

    if islogical(x_0) % Add the initial guess
        x = Q_previter(:,1:k-1) * y; % Compute x as sum_j y_j q_j - the last "just added" column of Q shouldnt be used
    else
        x = x_0 + Q_previter(:,1:k-1) * y; % Compute x as sum_j y_j q_j - the last "just added" column of Q shouldnt be used            
    end

timings = [TimeP1,TimeP2,k-1];
end




