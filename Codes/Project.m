% Loading the matrix file
load ash219.mat

M = Problem.A;
[m, n] = size(M);

% Generate RHS b if not available
% Create it as a perturbation of the row sum of the matrix
try
    b = Problem.b;
    disp('We use the given RHS vector.')
catch ME
   if (strcmp(ME.identifier,'MATLAB:nonExistentField'))
        disp('We create a RHS vector.')
        row = sum(M, 2);
        rand = - 5 + (7+5)*rand(m, 1);
        b = row + rand;
   end
end

% Step 1 - Compute LU factorization
% For stability, we use threshold pivoting
% For sparsity, we use the column reordering algorithm
p = colamd(M);
[L, U, P, Q] = lu(M(:, p), 0.1);
err = P*M*Q - L*U;
fprintf('The 1-norm of the difference between P*M*Q and L*U is %.3f.\n', ...
    norm(err, 1));

% Plot the sparsity patterns of matrix M, and L and U
figure
subplot(1,3,1)
spy(M)
title("Original Matrix")
subplot(1,3,2)
spy(L)
title('L factor')
subplot(1,3,3)
spy(U)
title('U factor')

% Step 2a - Split L and b into two parts
L1 = L(1:n, :);
L2 = L(n+1:end, :);

newb = P*b;
b1 = newb(1:n);
b2 = newb(n+1:end);

% Step 2b - Solve for c in L_1.c = b_1
c = L1\b1;

% Step 3 - Solve UQ^T x_1 = c;
z = U\c;
x1 = Q.'\z;  % dim(x1) = n by 1

% Step 4 - Form vector d
d = b2 - L2*c;  % dim(d) = (m-n) by 1
fprintf('The norm of the residual vector is %.3f.\n', norm(d, 2))

% Step 5 - Is residual norm below the (user-defined) threshold?
if norm(d, 2) < 10
    disp('Process stopped.')
    fprintf('Solution vector yields a residual of %.3f.\n', ...
        norm(b - M*x1, 2))
else
    disp('We will continue the process.')
    M1 = L.'*L;  % dim(sym_mat) = n by n

    % Find the Cholesky factorization of M1
    % Use the amd algorithm for reordering
    % Obtain a lower triangular Cholesky factor
    [H, flag, P] = chol(M1, 'lower');

    % Plot the sparsity patterns of matrix M1 and H
    figure
    subplot(1,2,1)
    spy(M1)
    title("Sparsity pattern of L'*L")
    subplot(1,2,2)
    spy(H)
    title('Sparsity pattern of chol(L''*L)')

    % Observe that HH.' = P.' * M1 * P
    % Hence, M1 = (PH)(PH).'
    % Now solve M1*z = L.'[0; d]
    if ~flag
        disp('Cholesky factorization successful.')
        % Form the RHS vector
        top = zeros(n, 1);
        rhs = L.'*[top; d];  % dim(rhs) = n by 1

        % Solve using Cholesky factors
        new_H = P*H;
        w = new_H\rhs;
        z = new_H.'\w;
 
        % Solve UQ^T x_2 = z
        temp = U\z;
        x2 = Q.'\temp;  % dim(x2) = n by 1

        % Find residual norm on final solution vector
        x = x1 + x2;
        disp('Process stopped.')
        fprintf('Solution vector yields a residual of %.3f.\n', ...
            norm(b - M*x, 2))
    else
        disp('Cholesky factorization failed.\n')
    end
end
