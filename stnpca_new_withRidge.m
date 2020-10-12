%-------------------------------Function Information-------------------------
% This function is used for extrcting Principle Components by Supervised 
% Tensor PCA with the ridge penelty
% Author: Wenlin Wu
% Created 2018-11-25


function[V U obj B beta alpha] = stnpca_n2(X, Y, k, S, gamma, delta)


%[V U lambda B alpha] = stnpca(X, Y, k, options)
%function used for extracting factors by supervised tensor PCA
%----------------------
%input information:
% X: PxPxN tensor, stack the PxP network matrix of N subjects
% Y: Nx1 vector, the response variables of N subjects
% k: int, the number of factors to extrct
% S: int, number of labeled data
% gamma: double, supervised information parameter
% delta: double, penalty information parameter
% v0, optinal, Px1 vector, the starting value of v
% u0, optional, Nx1 vector, the starting value of u
% maxit, optional, int, the maximum iteration, default:1000
% thr, optional, double, the stopping creteria, default:0.00001
%----------------------
%output information:
%V: Pxk matrix, factor matrix for Network mode
%U: Nxk matrix, factor matrix for Subject mode
%B: kx1 vector, regression coefficient
%alpha Nx1 vector, intercept
%beta coeff for u
%obj store value for objective function

N = size(X,3);
P = size(X,1);

% %input check
% p = inputParser;
% p.addRequired('X', @(x)validateattributes(x,{'tensor'},{}));
% p.addRequired('Y', @(y)validateattributes(y,{'double'},{'numel',N}));
% p.addRequired('k', @(y)validateattributes(y,{'int'},{}));
% 
% defaultu0 = randn(N,5); %start with 5 random point if not given start value
% defaultv0 = randn(P,5); %start with 5 random point if not given start value
% defaultmaxit = 1e3;
% defaultthr = 1e-5;
% 
% p.addOptional('v0', defaultv0);
% p.addOptional('u0', defaultu0);
% p.addOptional('maxit', defaultmaxit);
% p.addOptional('thr', defaultthr);
% 
% parse(X,Y,k,varargin{:});
% 
% X = p.Results.X;
% Y = p.Results.Y;
% k = p.Results.k;
% u0 = p.Results.u0;
% v0 = p.Results.v0;
% maxit = p.Results.maxit;
% thr = p.Results.thr;
% 
% size_u0 = size(u0,2);
% size_v0 = size(v0,2);

D = [eye(S),zeros(S,N-S)];
Rep = 5; %repeat and have the one with smallest objective function
maxit = 3000;
thr = 1e-7;
U = [];
V = [];
obj = [];
l = [];

for i = 1:k
    %save results for each iteration
    obj_all = [];
    v_store = [];
    u_store = [];
    ttt_store = [];
    if i==1
        for r1 = 1:Rep %local search find a best obj

            %initialize u,v
            u = randn(N, 1);
            u = u/norm(u);
            v = randn(P, 1);
            v = v/norm(v);
            alpha = randn(1);

            dif = 1; %assert whether the algorithm converge
            iter = 1;
            obj_store = [];
            u_store = [];
            while dif>thr && iter<maxit
                beta = (gamma/(gamma*(D*u)'*(D*u)+delta))*(D*u)'*(Y-alpha*ones(S,1));
                alpha = (1/S)*sum(Y-beta*(D*u));
                lambda = ttv(X,{v,v,u},[1,2,3]);
                %u = ((lambda*double(ttv(X,{v,v},[1,2])))+(beta*...
                   %D'*(Y-alpha*ones(S,1))));
                %u = u/norm(u);
%                 [evector,~] = eig(double(ttv(X,u,3)));
%                 v = real(evector(:,end));
                u_1 = ((lambda^2+delta*beta^2)*eye(size(D'*D,2)) + gamma*beta^2*(D'*D));
                u_2_1 = u_1\(lambda*double(ttv(X,{v,v},[1,2])));
                u_2_2 = u_1\(gamma*beta*D'*(Y-alpha*ones(S,1)));
                u = u_2_1+u_2_2;
                u_store = [u_store u_2_1 u_2_2];
                [v,~] = eigs(double(ttv(X,u,3)),1);
                obj1 = norm(X-full(ktensor(lambda,v,v,u)))^2;
                obj2 = gamma*norm(Y-(D*u)*beta-alpha*ones(S,1))^2;
                obj3 = delta*(norm(beta)^2);
                obj = obj1+obj2+obj3;
                obj_store = [obj_store obj];
                if iter == 1
                    dif = obj_store(end);
                else
                    dif = abs((obj_store(end)-obj_store(end-1))/obj_store(end-1));
                end
            iter = iter +1;
            end
            obj_all = [obj_all obj_store(end)];
            v_store = [v_store v];
            u_store = [u_store u];
        end
        [~,idx] = min(obj_all);
        v = v_store(:,idx);
        u = u_store(:,idx);
        V = [V v];
        U = [U u];
        obj = min(obj_all);
        lambda = ttv(X,{v,v,u},[1,2,3]);
        Xhat = X-full(ktensor(lambda,v,v,u));
        
 
    else
        %proj = eye(P) - V*V';
        for r = 1:Rep %local search find a best obj
            proj = eye(P) - V*V';
            %initialize u,v
            u = randn(N, 1);
            u = u/norm(u);
            v = randn(P, 1);
            v = v/norm(v);
            alpha = randn(1);
            B = randn(i-1,1);
            
            dif = 1; %assert whether the algorithm converge
            iter = 1;
            obj_store = [];
            while dif>thr && iter<maxit                
                lambda = ttv(Xhat,{proj*v,proj*v,u},[1,2,3]);
                %lambda = abs(ttv(Xhat,{proj*v,proj*v,u},[1,2,3]));
                beta = (gamma/(gamma*(D*u)'*(D*u)+delta))*(D*u)'*(Y-alpha*ones(S,1)-D*U*B);
                alpha = (1/S)*(sum((Y-beta*(D*u)-D*U*B)));
                B1 = gamma*(gamma*(D*U)'*(D*U)+delta*eye(size((D*U)'*(D*U),2)));
                B2 = (D*U)'*(Y-D*u*beta-alpha*ones(S,1));
                B = B1\B2;
%                 u = ((lambda*double(ttv(Xhat,{proj*v,proj*v},[1,2])))+(beta*...
%                    D'*(Y-alpha*ones(S,1)-D*U*B)));
%                 u = u/norm(u);
%                 [evector,~] = eig((proj)'*double(ttv(Xhat,u,3))*proj);
%                 v = real(evector(:,end));
                u_1 = ((lambda^2+delta*beta^2)*eye(size(D'*D,2)) + gamma*beta^2*(D'*D));
                u_2_1 = u_1\(lambda*double(ttv(Xhat,{proj*v,proj*v},[1,2])));
                u_2_2 = u_1\(gamma*beta*D'*(Y-alpha*ones(S,1)-D*U*B));
                u = u_2_1+u_2_2;
                [v,~] = eigs(proj*double(ttv(Xhat,u,3))*proj,1);
                obj1 = norm(Xhat-full(ktensor(lambda,proj*v,proj*v,u)))^2;
                obj2 = gamma*norm(Y-D*u*beta-D*U*B-alpha*ones(S,1))^2;
                obj3 = delta*(norm(beta)^2+norm(B)^2);
                obj = obj1+obj2+obj3;
                obj_store = [obj_store obj];
                l = [l lambda];
                if iter == 1
                    dif = obj_store(end);
                else
                    dif = abs((obj_store(end)-obj_store(end-1))/obj_store(end-1));
                end
                iter = iter +1;
            end
            obj_all = [obj_all obj_store(end)];
            v_store = [v_store v];
            u_store = [u_store u];
            plot(obj_store)
        end
        [~,idx] = min(obj_all);
        v = v_store(:,idx);
        u = u_store(:,idx);
        V = [V v];
        U = [U u];
        obj = min(obj_all);
        lambda = ttv(Xhat,{v,v,u},[1,2,3]); 
        %lambda = abs(ttv(Xhat,{v,v,u},[1,2,3]));
        Xhat = Xhat-full(ktensor(lambda,v,v,u));
        disp(i)
    end
end

end

    



