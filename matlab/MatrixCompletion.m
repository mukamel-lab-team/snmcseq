function x = MatrixCompletion(A, Mask, gamma, k)
% Matrix Completion using nuclear norm regularization and FISTA algorithm.
% Solves the optimization (1/2)(||Mask.*(x-A)||_2)^2+gamma*||x||_{*}.
% A is the incomplete matrix, Mask indicates the valid elements, x is the completed matrix using this algorithm.
fprintf('\n\n\n*******Matrix Completion with Nuclear Norm Regularization*******\n\n');

% Initialization.
xk = zeros(size(A)); vk = xk;
tk = 1; % Initial step sizes.
r1 = 0.8; r2 = 0.1;

maxIter = 5; count = 1;
Maincost = zeros(maxIter,1);
DataFidelity = zeros(maxIter,1);
Regularization = zeros(maxIter,1);

% Iterative loop
fprintf('Begin Optimization...\n\n');
tic

while (count<=maxIter)
    t = tk/r1;
    subiteration = 0;
    if count == 2
        r2 = 0.5;
    end
    while 1 % Lipschitz Gradient Condition
        if(count == 1)
            theta = 1;
        else
            r = tk/t;
            theta = thetak*(sqrt(4*r+thetak^2) - thetak)/(2*r);
        end
        
        y = (1-theta)*xk + theta*vk;
        delta_fy = (y-A).*Mask;
        
        x = proxNuclearNorm(y-t*delta_fy,t*gamma,k);
        
        % Compute f(y)
        fy = 1/2*sum(delta_fy(:).^2);
        UpperBound_x = fy + sum(sum(delta_fy.*(x-y))) + 1/(2*t)*sum(sum((x-y).^2));
        
        % Compute f(x)
        DataDiff_x = (x-A).*Mask;
        fx = 1/2*sum(DataDiff_x(:).^2);
        
        if fx <= UpperBound_x
            break
        end
        
        t = r2*t;
        subiteration = subiteration+1;
    end
    
    tk = t;
    thetak = theta;
    vk = xk + 1/theta*(x - xk);
    xk = x;
    
    DataFidelity(count) = fx;
    Regularization(count) = sum(svd(x));
    Maincost(count) = DataFidelity(count) + gamma*Regularization(count);
    
    if mod(count,5) == 0
        figure(3);semilogy(Maincost,'r');title(['Maincost iteration: ' num2str(count) ', gamma: ' num2str(gamma)]); pause(0.01);
        figure(4);semilogy(DataFidelity,'r');title(['Regularization iteration: ' num2str(count) ', gamma: ' num2str(gamma)]); pause(0.01);
        figure(5);semilogy(Regularization,'r');title(['Data Fidelity iteration: ' num2str(count) ', gamma: ' num2str(gamma)]); pause(0.01);
    end
    count = count+1;
end

fprintf('\n\nEnd Optimization.');
fprintf('\n\n\n*******Matrix Completion with Nuclear Norm Regularization*******\n\n');

end



function prox = proxNuclearNorm(x,t,k)
% Solves the proximal operator of Nuclear Norm
% Solves the optimization (1/(2t))(||x-y||_2)^2 + ||x||_{*}.
% Keep k largest singular values then perform thresholding

[u,s,v]=svds(x,k);
diags = diag(s)-t;
s = diag(max(diags,0));
prox = u*s*v';

end
      

