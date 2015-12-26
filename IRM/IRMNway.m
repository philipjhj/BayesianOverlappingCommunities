function [L,cpu_time,Z,eta,sample,West]=IRMNway(A,W,noc,opts)

% Non-parametric IRM of N-partite graphs based on collapsed Gibbs sampling with split
% merge moves. The method includes sampling of hyper-parameters.
%
% Usage: 
%  [L,cpu_time,Z1,Z2,eta,sample,West]=IRMBipartite(A,W,noc,opts)
%
% Input:
%   A       N-way array
%   W       N-way binary matrix of missing values
%   noc     1 x N-modes vector indicating the number of initial clusters in Z{n}
%   opts.
%           maxiter     maximum number of iterations
%           Z{n}        noc x I(n) matrix of community identities
%           dZstep      number of iteratiosn between each recorded sample
%           verbose     1: display iteration results, 0: no display
%           eta0        1x2 vector of intial values of pseudo link and non-link counts between the groups (default: [1 1])
%           type        'Binary' (default), 'Weighted'
% Output:
%   L           Log-likelihood function at each iteration
%   cpu_time    cpu-time cost for each iteration
%   Z           cell array of estimated assigment matrices
%   sample      sampled parameters at each dZstep iterationCluster
%   eta         Average between group link densitites
%   West        Estimated link values for missing links and non-links
%               for type='Categorical' this corresponds to the probability
%               of generating the observed class.
% 
% This code is provided as is without any kind of warranty!
%
% Written by Morten Mørup
 
if nargin<4
    opts=struct;
end

% Make sure links treated as missing are removed from estimation
par.type=mgetopt(opts,'type','Binary');
N=size(A);
nmodes=length(N);
if strcmp(par.type,'Binary')
    A=logical(A-A.*W); % Remove missing links
else
     A=A-A.*W; % Remove missing links
end
nnzW=nnz(W);

% Initialize Z
if isfield(opts,'Z')
    Z=opts.Z;
else    
    for n=1:ndims(A)
        ind=ceil(noc(1)*rand(1,N(n)));
        Z{n}=sparse(ind,1:N(n),ones(1,N(n)),noc(1),N(n));
        Z{n}=mgetopt(opts,'Z',Z{n});
        Z{n}=full(Z{n});
        Z{n}(sum(Z{n},2)==0,:)=[];
    end
end

% Set remaining parameters
if isfield(opts,'alpha')
    alpha=opts.alpha;
else
    alpha=log(size(A));
end
maxiter=mgetopt(opts,'maxiter',50);
verbose=mgetopt(opts,'verbose',1);
switch par.type
    case 'Binary'
        eta0=mgetopt(opts,'eta0',[1 1]); % pseudo counts of links and non-links between clusters
    case 'Weighted'
        eta0=mgetopt(opts,'eta0',[1 1e-6]); % pseudo counts of links and non-links between clusters
end
sample_step=mgetopt(opts,'sample_step',25); % Steps between samples
nsampleiter=mgetopt(opts,'nsampleiter',25);

% Set algorithm variables
sample=struct;
L=zeros(1,maxiter);
cpu_time=L;

West=0;
sstep=0;
westiter=0;
Q=-inf;
Qbest=-inf;
iter=0;

str2b=[];
str2a=['%' num2str(nmodes*12-5*nmodes) 's |'];
for k=1:nmodes
    str(((k-1)*8)+1:((k-1)*8)+8)='%12.0f |';        
    str2b=[str2b ['noc' num2str(k) ' ']];
end
if verbose % Display algorithm    
    disp(['Non-parametric N-wa clustering based on the IRM model for Binary or Weighted graphs'])
    dheader = sprintf(['%12s | %12s | %12s | ' str2a ' %12s | %12s'],'Iteration','logP','dlogP/|logP|',str2b,'time');
    dline = sprintf('-------------+--------------+--------------+--------------+--------------+--------------');
    disp(dline);
    disp(dheader);
    disp(dline);
end


% Main Loop
while iter<maxiter
   
    iter=iter+1;
    tic;
    Qold=Q;
    
    
    % Gibbs sampling of Z for each mode 
    for k=1:ndims(A) 
        modes=setdiff(1:ndims(A),k);
        Zm=1;
        for kk=modes
            Zm=kron(Z{kk},Zm);
        end        
        Am=reshape(permute(A,[k, modes]),[N(k),prod(N(modes))]);            
        Wm=reshape(permute(W,[k, modes]),[N(k),prod(N(modes))]);         
        par.ZA=Zm*Am';
        par.ZW=Zm*Wm';
    
        [Z{k}, logP_A, logP_Z(k)]=Gibbs_sample_ZIRM(Z{k},Zm,Am,Wm,eta0,alpha(k),randperm(N(k)),1,par);            
        noc(k)=size(Z{k},1);
        for t=1:10
            [Z{k},logP_A,logP_Z(k)]=split_merge_sample_Z(Z{k},Zm,Am,Wm,eta0,alpha(k),logP_A,logP_Z(k),1,par);
        end  
        [~,ind]=sort(sum(Z{k},2),'descend');
        Z{k}=Z{k}(ind,:);
        noc(k)=size(Z{k},1);
        
        % sample alpha        
        [logP_Z(k),alpha(k)]=sample_alpha(sum(Z{k},2),N(k),alpha(k));        
    end    
    
    % sample eta0        
    [logP_A,eta0]=sample_eta0(Z{k},Zm,eta0,par,Am,Wm);           
    etam=calculateEta(Zm,Z{k},eta0,par);
    eta=permute(reshape(etam,[noc(k), noc(modes)]),[2:nmodes 1]);    
    
    % Evaluate result
    Q=logP_A+sum(logP_Z);
    dQ=Q-Qold;    
    L(iter)=Q;
    t_iter=toc;
    cpu_time(iter)=t_iter;    
   
    % Display iteration
    if rem(iter,1)==0 && verbose        
        disp(sprintf(['%12.0f | %12.4e | %12.4e |' str ' %12.4f '],iter,Q,dQ/abs(Q),noc,t_iter));
    end
    
    % Ztore sample
    if mod(iter,sample_step)==0         
        sstep=sstep+1;
        sample.iteration(sstep)=iter;
        sample.Z{sstep}=Z;                
        sample.eta{sstep}=eta;
        sample.alpha{sstep}=alpha;
        sample.eta0{sstep}=eta0;
    end        
    if Q>Qbest
       sample.MAP.Z=Z;        
       sample.MAP.iter=iter;
       sample.MAP.L=Q;
       sample.MAP.eta=eta;
       sample.MAP.eta0=eta0;
       sample.MAP.alpha=alpha;
       Qbest=Q;
    end
        
        
    % Estimate missing link probability of sample
    if iter>maxiter-nsampleiter && nnzW>0
        westiter=westiter+1;        
        if iter==maxiter-nsampleiter+1
            disp(['Initiating estimation of missing links for the last ' num2str(nsampleiter) ' iterations']);   
        end      
        step=10000;
        [Iw,Jw]=find(Wm);
        val=zeros(1,length(Iw));            
        for kk=1:ceil((length(Iw)/step))
              ind=(kk-1)*step+1:min([kk*step, length(Iw)]);   
              val(ind)=sum(Z{k}(:,Iw(ind)).*(etam(:,:)*Zm(:,Jw(ind))))+eps;
        end                              
        West=West+permute(reshape(full(sparse(Iw,Jw,val,size(Am,1),size(Am,2))),[N(k), N(modes)]),[2:nmodes 1]);      
    end    

end

% Average link predictions
if nnzW>0
    West=West/westiter;    
end
% Display final iteration
if verbose   
  disp('Result of final iteration');
  disp(sprintf('%12.0f | %12.4e | %12.4e | %12.0f | %12.0f |%12.4f ',iter,Q,dQ/abs(Q),noc(1),noc(2),t_iter));
end

% -------------------------------------------------------------------------
% Parser for optional arguments
function var = mgetopt(opts, varname, default, varargin)
if isfield(opts, varname)
    var = getfield(opts, varname); 
else
    var = default;
end
for narg = 1:2:nargin-4
    cmd = varargin{narg};
    arg = varargin{narg+1};
    switch cmd
        case 'instrset',
            if ~any(strcmp(arg, var))
                fprintf(['Wrong argument %s = ''%s'' - ', ...
                    'Using default : %s = ''%s''\n'], ...
                    varname, var, varname, default);
                var = default;
            end
        otherwise,
            error('Wrong option: %s.', cmd);
    end
end

%-------------------------------------------------------------------------
function eta=calculateEta(Z2,Z1,eta0,par)
        noc1=size(Z1,1);
        noc2=size(Z2,1);
        sumZ1=sum(Z1,2);
        sumZ2=sum(Z2,2);        
        eta=zeros(noc1,noc2);
        ZAZt=Z1*par.ZA';
        ZWZt=Z1*par.ZW';            
        switch par.type
            case 'Binary'
                ZZ=sumZ1*sumZ2'-ZAZt-ZWZt;                
                n_link=ZAZt+eta0(1);                        
                n_nonlink=ZZ+eta0(2);   
                eta=n_link./(n_link+n_nonlink);                                
            case 'Weighted'                
                ZZ=sumZ1*sumZ2'-ZWZt;                
                n_link=ZAZt+eta0(1);                        
                n_nonlink=ZZ+eta0(2);   
                eta=n_link./(n_nonlink);        
        end
        

% -------------------------------------------------------------------------  
function [Z2,logP_A,logP_Z2]=split_merge_sample_Z(Z2,Z1,A,W,eta0,alpha,logP_A,logP_Z2,mode,par)
    %[logP_A_t,logP_Z2_t]=evalLikelihood(Z2,Z1,A,W,eta0,alpha,mode);                                 

    noc2=size(Z2,1);
    if mode==2
        J=size(A,2);    
    else
        J=size(A,1);    
    end
    
    % step 1 select two observations i and j        
    ind1=ceil(J*rand);        
    ind2=ceil((J-1)*rand);
    if ind1<=ind2
       ind2=ind2+1;
    end
    clust1=find(Z2(:,ind1));
    clust2=find(Z2(:,ind2));

    if clust1==clust2 % Split   
        setZ=find(sum(Z2([clust1 clust2],:)));    
        setZ=setdiff(setZ,[ind1,ind2]);
        n_setZ=length(setZ);
        Z2_t=Z2;
        Z2_t(clust1,:)=0;        
        comp=[clust1 noc2+1];               
        Z2_t(comp(1),ind1)=1;
        Z2_t(comp(2),ind2)=1;

        % Reassign by restricted gibbs sampling        
        if n_setZ>0
            for rep=1:3
                [Z2_t,logP_A_t,logP_Z2_t,logQ_trans,comp]=Gibbs_sample_ZIRM(Z2_t,Z1,A,W,eta0,alpha,setZ(randperm(n_setZ)),mode,par,comp);                     
            end     
        else
           logQ_trans=0;
           [logP_A_t,logP_Z2_t]=evalLikelihood(Z2_t,Z1,A,W,eta0,alpha,mode,par);                 
        end
        
        % Calculate Metropolis-Hastings ratio
        a_split=rand<exp(logP_A_t+logP_Z2_t-logP_A-logP_Z2-logQ_trans);         
        if a_split
           disp(['Splitting cluster'])
           logP_A=logP_A_t;
           logP_Z2=logP_Z2_t;
           Z2=Z2_t;
        end
    else % Merge                                     
        Z2_t=Z2;
        Z2_t(clust1,:)=Z2_t(clust1,:)+Z2_t(clust2,:);
        setZ=find(Z2_t(clust1,:));           
        Z2_t(clust2,:)=[];        
        if clust2<clust1
            clust1_t=clust1-1;
        else 
            clust1_t=clust1;
        end
        noc2_t=noc2-1;

        % calculate likelihood of merged cluster       
        [logP_A_t,logP_Z2_t]=evalLikelihood(Z2_t,Z1,A,W,eta0,alpha,mode,par);                

        % Zplit the merged cluster and calculate transition probabilties                        
        setZ=setdiff(setZ,[ind1,ind2]);
        n_setZ=length(setZ);
        Z2_tt=Z2_t;
        Z2_tt(clust1_t,:)=0;        
        comp=[clust1_t noc2_t+1];               
        Z2_tt(comp(1),ind1)=1;
        Z2_tt(comp(2),ind2)=1;                
        
        % Reassign by restricted gibbs sampling
        if n_setZ>0
            for rep=1:2        
                [Z2_tt,logP_A_tt,logP_Z2_tt,logQ_trans,comp]=Gibbs_sample_ZIRM(Z2_tt,Z1,A,W,eta0,alpha,setZ(randperm(n_setZ)),mode,par,comp);                
            end
            Force=[1 2]*Z2([clust1 clust2],:);        
            [Z2_tt,logP_A_tt,logP_Z2_tt,logQ_trans]=Gibbs_sample_ZIRM(Z2_tt,Z1,A,W,eta0,alpha,setZ(randperm(n_setZ)),mode,par,comp,Force);                        
        else
            logQ_trans=0;                   
        end
        a_merge=rand<exp(logP_A_t+logP_Z2_t-logP_A-logP_Z2+logQ_trans); 
        
        if a_merge
          disp(['Merging cluster'])
          logP_A=logP_A_t;
          logP_Z2=logP_Z2_t;
          Z2=Z2_t;          
        end
    end

% -------------------------------------------------------------------------  
function [logP_A,logP_Z2,ZAZt,ZWZt,sumZ2]=evalLikelihood(Z2,Z1,A,W,eta0,alpha,mode,par)

    switch par.type
        case 'Binary'
            my_func=@betaln;
            const=my_func(eta0(1),eta0(2));   
        case 'Weighted'
            my_func=@poissonln;
            const=my_func(eta0(1),eta0(2));           
    end
    if mode==2
        J=size(A,2);
    else
        J=size(A,1);
    end
    noc1=size(Z1,1);
    noc2=size(Z2,1);              
    
    if mode==2
        ZAZt=Z1*A*Z2';
        ZWZt=Z1*W*Z2';
    else
        ZAZt=Z1*A'*Z2';
        ZWZt=Z1*W'*Z2';
    end    

    sumZ1=sum(Z1,2);
    sumZ2=sum(Z2,2);
    switch par.type
        case 'Binary'
            ZZ=sumZ1*sumZ2'-ZAZt-ZWZt;                
        case 'Weighted'
            ZZ=sumZ1*sumZ2'-ZWZt;                        
    end
    n_link=ZAZt+eta0(1);    
    n_nonlink=ZZ+eta0(2);         
    logP_A=sum(sum(my_func(n_link,n_nonlink)))-noc1*noc2*const;
    logP_Z2=noc2*log(alpha)+sum(gammaln(full(sumZ2)))-gammaln(J+alpha)+gammaln(alpha);      
    
% -------------------------------------------------------------------------
function [Z2,logP_A,logP_Z2,logQ_trans,comp]=Gibbs_sample_ZIRM(Z2,Z1,A,W,eta0,alpha,JJ,mode,par,comp,Force)    
    
    if nargin<11
        Force=[];
    end
    if nargin<10
        comp=[];
    end    
    logQ_trans=0;

    switch par.type
        case 'Binary'
            my_func=@betaln; 
            const=my_func(eta0(1),eta0(2));                   
        case 'Weighted' 
            my_func=@poissonln;
            const=my_func(eta0(1),eta0(2));                           
    end
    if mode==2
        [I,J]=size(A);
    else
        [J,I]=size(A);
    end      
        
    t=0;   
    sumZ1=sum(Z1,2);
    noc1=length(sumZ1);    
    sumZ2=sum(Z2,2);
    noc2=length(sumZ2);        
    
    ZA=zeros(noc1,J);
    ZW=zeros(noc1,J);
    n_link=zeros(noc1,noc2);
    n_nonlink=n_link;
    
    ZA=par.ZA;
    ZW=par.ZW;
    
    ZAZt=ZA*Z2';
    ZWZt=ZW*Z2';                    
    n_link(:,:)=ZAZt+eta0(1);
    switch par.type
        case 'Binary'
            n_nonlink(:,:)=sumZ1*sumZ2'-ZAZt-ZWZt+eta0(2);        
        case 'Weighted'
            n_nonlink(:,:)=sumZ1*sumZ2'-ZWZt+eta0(2);        
    end
    
    beta_eval=my_func(n_link,n_nonlink);    
    
    for k=JJ           
        t=t+1;
        if mod(t,5000)==0
            disp(['sampling ' num2str(t) ' out of ' num2str(J) ' nodes']);
        end
        
        % Remove effect of Z2(:,k)        
        sumZ2=sumZ2-Z2(:,k);   
        Z1AZ2k=ZA(:,k);
        Z1WZ2k=ZW(:,k);  
        switch par.type
            case 'Binary'
                nZ1AZ2k=sumZ1-Z1AZ2k-Z1WZ2k;
            case 'Weighted'
                nZ1AZ2k=sumZ1-Z1WZ2k;
        end
        d=find(Z2(:,k));        
        if ~isempty(d)
            n_link(:,d)=n_link(:,d)-Z1AZ2k;                                            
            n_nonlink(:,d)=n_nonlink(:,d)-nZ1AZ2k;                                                       
            Z2(:,k)=0;               
        end
        
        if isempty(comp)
            if sumZ2(d)==0 % Remove singleton cluster            
                v=1:noc2;
                v(d)=[];
                d=[];
                noc2=noc2-1;                               
                P=sparse(1:noc2,v,ones(1,noc2),noc2,noc2+1);                        
                Z2=P*Z2;                        
                sumZ2=sumZ2(v,1);            
                n_link=n_link(:,v);            
                n_nonlink=n_nonlink(:,v);            
                beta_eval=beta_eval(:,v);            
            end                                       

            % Calculate probability for existing communties as well as proposal cluster                                                                 
            beta_eval(:,d)=my_func(n_link(:,d),n_nonlink(:,d)); % removed the constant -my_func(Ap,An)))                                           
            sum_beta_eval=sum(sum(beta_eval));
            e=ones(1,noc2);
            sum_beta_eval_d=sum_beta_eval-[sum(beta_eval,1) 0];
            TTT1=cat(2,n_link+Z1AZ2k(:,e),Z1AZ2k+eta0(1));
            TTT2=cat(2,n_nonlink+nZ1AZ2k(:,e),nZ1AZ2k+eta0(2));
            beta_eval_d=my_func(TTT1,TTT2);                                        
            logQ=sum_beta_eval_d'+[sum(beta_eval_d(:,1:noc2,:),1) (sum(beta_eval_d(:,noc2+1),1)-noc1*const)]'; % removed the constant -my_func(Ap,An)))                                             

                        
            % Zample from posterior                        
            QQ=exp(logQ-max(logQ));
            weight=[sumZ2; alpha];
%             % code test
%             noc1t=ceil(size(Z2,1)*rand);
%             noc2t=ceil((size(Z2,1)+1)*rand);
%             Z2_t=Z2;
%             Z2_t(noc1t,k)=1;
%             [logP_A_t1,logP_Z2_t1]=evalLikelihood(Z2_t,Z1,A,W,eta0,alpha,mode,par);     
%             Z2_t=Z2;
%             Z2_t(noc2t,k)=1;
%             [logP_A_t2,logP_Z2_t2]=evalLikelihood(Z2_t,Z1,A,W,eta0,alpha,mode,par);                 
%             a1=logP_A_t1+logP_Z2_t1-(logP_A_t2+logP_Z2_t2);                        
%             a2=logQ(noc1t)+log(weight(noc1t))-(logQ(noc2t)+log(weight(noc2t)));
%             difference=a1-a2;
%             if abs(difference)>1e-9
%                 keyboard;
%             end
%              
            
            
            QQ=weight.*QQ;        
            ind=find(rand<cumsum(QQ/sum(QQ)),1,'first');     
            Z2(ind,k)=1;   
            if ind>noc2            
                noc2=noc2+1;
                sumZ2(noc2,1)=0;
                n_link(:,noc2)=eta0(1);                                
                n_nonlink(:,noc2)=eta0(2);                     
                beta_eval(:,noc2)=0;                    
            end
            beta_eval_d=beta_eval_d(:,ind);            
        else            
            % Calculate probability for existing communties as well as proposal cluster                                                                  
            beta_eval(:,d)=my_func(n_link(:,d),n_nonlink(:,d)); % removed the constant -my_func(Ap,An)))                                           
            sum_beta_eval=sum(sum(beta_eval));                        
            e=ones(1,2);
            sum_beta_eval_d=sum_beta_eval-sum(beta_eval(:,comp),1);
            beta_eval_d=my_func(n_link(:,comp)+Z1AZ2k(:,e),n_nonlink(:,comp)+nZ1AZ2k(:,e));                                        
            logQ=sum_beta_eval_d'+sum(beta_eval_d,1)'; % removed the constant -my_func(Ap,An)))                                             
           
            % Zample from posterior                        
            QQ=exp(logQ-max(logQ));
            weight=sumZ2(comp);
%                         % code test
%             noc1t=comp(1);
%             noc2t=comp(2);
%             Z2_t=Z2;
%             Z2_t(noc1t,k)=1;
%             [logP_A_t1,logP_Z2_t1]=evalLikelihood(Z2_t,Z1,A,W,eta0,alpha,mode,par);     
%             Z2_t=Z2;
%             Z2_t(noc2t,k)=1;
%             [logP_A_t2,logP_Z2_t2]=evalLikelihood(Z2_t,Z1,A,W,eta0,alpha,mode,par);                 
%             a1=logP_A_t1+logP_Z2_t1-(logP_A_t2+logP_Z2_t2);                        
%             a2=logQ(1)+log(weight(1))-(logQ(2)+log(weight(2)));
%             difference=a1-a2;
%             if abs(difference)>1e-9
%                 keyboard;
%             end
%            
            
            QQ=weight.*QQ;
            QQ=QQ/sum(QQ);            
            if isempty(Force)
                ind=find(rand<cumsum(QQ),1,'first');
            else 
                ind=Force(k);
            end
            q_tmp=logQ-max(logQ)+log(weight);
            q_tmp=q_tmp-log(sum(exp(q_tmp)));            
            logQ_trans=logQ_trans+q_tmp(ind);
            Z2(comp(ind),k)=1;  
            beta_eval_d=beta_eval_d(:,ind);            
            ind=comp(ind);
        end
                        
        % Re-enter effect of new s_k        
        sumZ2=sumZ2+Z2(:,k);
        n_link(:,ind)=n_link(:,ind)+Z1AZ2k;                                
        n_nonlink(:,ind)=n_nonlink(:,ind)+nZ1AZ2k;        
        beta_eval(:,ind)=beta_eval_d;
                
        % Remove empty clusters        
        if ~all(sumZ2)
            d=find(sumZ2==0);
            if ~isempty(comp)
                ind_d=find(d<comp);
                comp(ind_d)=comp(ind_d)-1;
            end
            v=1:noc2;
            v(d)=[];
            noc2=noc2-length(d);                               
            P=sparse(1:noc2,v,ones(1,noc2),noc2,noc2+1);                        
            Z2=P*Z2;                        
            sumZ2=sumZ2(v,1);            
            n_link=n_link(:,v);
            n_nonlink=n_nonlink(:,v);
            beta_eval=beta_eval(:,v);            
        end              
    end      
    % Calculate Likelihood for sampled solution     
    logP_Z2=noc2*log(alpha)+sum(gammaln(full(sumZ2)))-gammaln(J+alpha)+gammaln(alpha);
    logP_A=sum(sum(beta_eval))-noc1*noc2*const;     
    
    % ---------------------
    function C=poissonln(A,B)
        C = gammaln(A)-A.*log(B);
        
    %--------------------------------------------------------------------
function [logZ,alpha]=sample_alpha(sumZ,N,alpha)

max_iter=100;
K=length(sumZ);
const=sum(gammaln(sumZ));
logZ=K*log(alpha)+const-gammaln(N+alpha)+gammaln(alpha);
accept=0;
for sample_iter=1:max_iter  
    alpha_new=exp(log(alpha)+0.1*randn); % symmetric proposal distribution in log-domain (use change of variable in acceptance rate alpha_new/alpha) 
    logZ_new=K*log(alpha_new)+const-gammaln(N+alpha_new)+gammaln(alpha_new); 
    if rand<(alpha_new/alpha*exp(logZ_new-logZ))
        alpha=alpha_new;
        logZ=logZ_new;
        accept=accept+1;       
    end
 end
 %disp(['accepted ' num2str(accept) ' out of ' num2str(max_iter) ' samples for alpha']);
    
 
 %--------------------------------------------------------------------
function [logP,eta0]=sample_eta0(Z1,Z2,eta0,par,A,W) 
       
max_iter=20;
noc1=size(Z1,1);
noc2=size(Z2,1);
sumZ1=sum(Z1,2);
sumZ2=sum(Z2,2);
ZAZt=Z1*par.ZA';
ZWZt=Z1*par.ZW';            
switch par.type    
    case 'Binary'
        my_func=@betaln; 
        const=my_func(eta0(1),eta0(2));                   
        ZZ=sumZ1*sumZ2'-ZAZt-ZWZt;                
        n_link_noeta0=ZAZt;                        
        n_nonlink_noeta0=ZZ;           
    case 'Weighted'
        my_func=@poissonln;
        const=my_func(eta0(1),eta0(2));                           
        ZZ=sumZ1*sumZ2'-ZWZt;                
        n_link_noeta0=ZAZt;                        
        n_nonlink_noeta0=ZZ+eta0(2);                 
end

n_link=n_link_noeta0+eta0(1);
n_nonlink=n_nonlink_noeta0+eta0(2);
logP=sum(sum(my_func(n_link,n_nonlink)))-noc1*noc2*const;     
accept=0;

for sample_iter=1:max_iter           
    eta0_new=exp(log(eta0)+0.1*randn(1,2)); % symmetric proposal distribution in log-domain (use change of variable in acceptance rate alpha_new/alpha)                         
    const_new=my_func(eta0_new(1),eta0_new(2));                                   
    n_link_new=n_link_noeta0+eta0_new(1);
    n_nonlink_new=n_nonlink_noeta0+eta0_new(2);
    logP_new=sum(sum(my_func(n_link_new,n_nonlink_new)))-noc1*noc2*const_new;     
    
%     % Test code
%     [logP_A_tmp,logP_Z2,ZAZt,ZWZt,sumZ2]=evalLikelihood(Z2,Z1,A,W,eta0,1,2,par);
%     [logP_A_tmp_new,logP_Z2,ZAZt,ZWZt,sumZ2]=evalLikelihood(Z2,Z1,A,W,eta0_new,1,2,par);
%     logP_A_tmp_new-logP_A_tmp
%     logP_new-logP
    
    if rand<(eta0_new(1)/eta0(1)*eta0_new(2)/eta0(2)*exp(logP_new-logP))
        eta0=eta0_new;
        logP=logP_new;        
        accept=accept+1;       
    end
end
 disp(['accepted ' num2str(accept) ' out of ' num2str(max_iter) ' samples for eta0']);
        