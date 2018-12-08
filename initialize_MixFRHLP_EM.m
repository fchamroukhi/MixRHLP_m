function param = initialize_MixFRHLP_EM(Y, G , K, Xbeta, XW, variance_type, init_kmeans, try_algo)

%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
%%%%%%%%%%%%%%%%%% Faicel Chamroukhi

[n, m]=size(Y);
p = size(Xbeta,2)-1;
q = size(XW,2)-1;


% % 1. Initialization of cluster weights
param.alpha_g=1/G*ones(G,1);
%2. Initialization of the model parameters for each cluster: W (pi_jgk), betak and sigmak
[W, pijgk] = init_hlp(G, K, q, XW, try_algo);
param.Wg = W;
param.pi_jgk = pijgk;
% betagk and sigmagk
if init_kmeans
    D = Y;
    max_iter_kmeans = 400;
    n_tries_kmeans = 20;
    verbose_kmeans = 0;
    
    res_kmeans = Kmeans_faicel(D,G,n_tries_kmeans, max_iter_kmeans, verbose_kmeans);
    
    for g=1:G
        Yg = D(res_kmeans.klas==g ,:); %if kmeans

        param_init =  init_regression_param(Yg, K,Xbeta, variance_type, try_algo);

        param.beta_g(:,:,g) = param_init.betak;
        if strcmp(variance_type,'common')
            param.sigma_g(g) = param_init.sigma;
        else
            param.sigma_g(:,g) = param_init.sigmak;
        end
    end
else
    ind = randperm(n);
    D=Y;
    for g=1:G
        if g<G
            Yg = D(ind((g-1)*round(n/G) +1 : g*round(n/G)),:);
        else
            Yg = D(ind((g-1)*round(n/G) +1 : end),:);
        end
        
        param_init =  init_regression_param(Yg, K, Xbeta, variance_type, try_algo);
        
        param.beta_g(:,:,g) = param_init.betak;
        
        if strcmp(variance_type,'common')
            param.sigma_g(g) = param_init.sigma;
        else
            param.sigma_g(:,g) = param_init.sigmak;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
function [Wg,  pi_jgk] =  init_hlp(G, K, q, XW, try_EM)

% init_hlp initialize the Hidden Logistic Process
%
%
%
%%%%%%%%%%%% FC %%%%%%%
% %   1. Initialisation de W (pi_jgk)

[nm, q1] = size(XW);

pi_jgk = zeros(nm,K,G);

Wg = zeros(q+1,K-1,G);
if  try_EM ==1
    for g=1:G
        %Wg(:,:,g) = zeros(q+1,K-1);%initialisation avec le vect null du vercteur param???tre du IRLS
        pi_jgk(:,:,g) = logit_model(Wg(:,:,g),XW);
    end
else
    for g=1:G
        Wg(:,:,g) = rand(q+1,K-1);%initialisation al??atoire du vercteur param???tre du IRLS
        pi_jgk(:,:,g) = logit_model(Wg(:,:,g),XW);
    end
end


%%%%%%%%%%%%%%%%%%%%%%
function para = init_regression_param(Y, K, phi,type_variance, try_EM)
% init_regression_param initialize the Regresssion model with Hidden Logistic Process
%
%
%
%%%%%%%%%%%%%%%%%%%% Faicel Chamroukhi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n, m] = size(Y);

if strcmp(type_variance,'common')
    s=0;
end
%

if try_EM ==1
    
    %decoupage de l'echantillon (signal) en K segments
    zi = round(m/K)-1;
    for k=1:K
        i = (k-1)*zi+1;
        j = k*zi;
        
        Yij = Y(:,i:j);
        Yij = reshape(Yij',[],1);

        X_ij = repmat(phi(i:j,:),n,1);

        bk = inv(X_ij'*X_ij)*X_ij'*Yij;
        
        para.betak(:,k) = bk;
        
        if strcmp(type_variance,'common')
            para.sigma = var(Yij);%1000;
        else
            mk = j-i+1 ;%length(Yij);

            z = Yij-X_ij*bk;
            sk = z'*z/(n*mk);
            para.sigmak(k) = sk;
            %para.sigmak(k) = var(Yij);
        end
    end
else % initialisation aleatoire
    Lmin= round(m/(K+1));%nbr pts min dans un segments
    tk_init = zeros(1,K+1);
    tk_init(1) = 0;
    K_1=K;
    for k = 2:K
        K_1 = K_1-1;
        temp = tk_init(k-1)+Lmin:m-K_1*Lmin;
        ind = randperm(length(temp));
        tk_init(k)= temp(ind(1));
    end
    tk_init(K+1) = m;
    
    %model.tk_init = tk_init;
    for k=1:K
        i = tk_init(k)+1;
        j = tk_init(k+1);
        Yij = Y(:,i:j);
        Yij = reshape(Yij',[],1);
        
        X_ij=repmat(phi(i:j,:),n,1);
        
        bk = inv(X_ij'*X_ij)*X_ij'*Yij;
        para.betak(:,k) = bk;
        
        if strcmp(type_variance,'common')
            para.sigma = var(Yij);%
        else
            mk = j-i+1 ;%length(Yij);
            z = Yij-X_ij*bk;
            sk = z'*z/(n*mk);
            para.sigmak(k) = sk;
            %para.sigmak(k) = var(Yij);
        end
        
    end
end
