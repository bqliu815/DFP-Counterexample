classdef PolarDFP
% SPDX-License-Identifier: MIT
% Exact prescribed polar angles and finite interpolation for classical DFP.
methods (Static)
function o=orbit(N,J,nuRatio)
    a=2*pi/9; nu=nuRatio*a;
    z=nu*(2*a-nu)/(2*nu-a); w=nu*(a+nu)/(2*nu-a);
    K=(2*nu-a)/(nu*(2*a-nu)*(a+nu)); L=4*a/(9*nu);
    X=[K,K*w;K*w,1+K*w*w];
    n=2*N; t0=J^(-1/3); S=diag([t0^2,1]);
    o.h0=S*(X\S); o.h0=(o.h0+o.h0')/2;
    o.J=J; o.nuRatio=nuRatio; o.cycles=N;
    o.targets=struct('A',K,'QE',z+2*a,'QO',w+a,'D0',2*a*(z+2*a)-2*a*a,'D1',a*(w+a)-a*a/2,'L',L);
    o.t=(J+(0:N)') .^(-1/3); o.A=zeros(N+1,1);o.A(1)=K;
    o.logW=zeros(n,1);o.delta=zeros(n,1);o.Q=zeros(n,1);
    o.armijo=zeros(n,1);o.strong=zeros(n,1);o.secant=zeros(n,1);o.update=zeros(n,1);
    o.zeta=ones(n+1,1);o.minEig=Inf;o.minSigma=Inf;
    for j=1:N
        t=o.t(j);t2=t*t;
        total=2*pi/t*expm1(log1p(t^3)/3);
        dz=expm1(-L*(K/X(1,1))^2*t^3);
        factors=[1,1+dz,1]; o.zeta(2*j)=1+dz;
        for p=1:2
            k=2*j-2+p;delta=total*(3-p)/3;
            zet=factors(p);zhat=factors(p+1);
            if p==1,change=dz;else,change=-dz;end
            Q=X(2,2)/X(1,2);q=t2*Q;
            vb=[-Q;1]; wb=[t2*change*cot(delta)-zet*t2*t2*Q;zhat];
            sigma=vb'*wb;
            P=eye(2)-wb*vb'/sigma;
            Y=P*X*P'+wb*wb'/sigma;Y=(Y+Y')/2;
            sd=sin(delta);cd=cos(delta);
            Rbar=[cd,t2*sd;-sd/t2,cd];
            Xnext=Rbar*Y*Rbar';Xnext=(Xnext+Xnext')/2;
            o.minEig=min(o.minEig,min(eig(Xnext)));o.minSigma=min(o.minSigma,sigma);
            assert(sigma>0 && Xnext(1,2)>0 && det(Xnext)>0,'Invalid polar state');
            v=[-q;1]; ww=[wb(1)/t2;wb(2)];
            H=diag([t2,1])*(X\diag([t2,1]));
            Hnew=diag([t2,1])*(Y\diag([t2,1]));
            Hy=H*ww;
            Hdfp=H-Hy*Hy'/(ww'*Hy)+v*v'/sigma;
            o.secant(k)=norm(Hnew*ww-v)/norm(v);
            o.update(k)=norm(Hdfp-Hnew,'fro')/norm(Hnew,'fro');
            wm1=-2*sin(delta/2)^2+q*sd;
            assert(wm1>0,'Radius did not decrease');
            W=1+wm1;
            o.logW(k)=log1p(wm1);o.delta(k)=delta;o.Q(k)=Q;
            o.armijo(k)=wm1*(2+wm1)/(2*zet*q*sd*W);
            o.strong(k)=zhat/zet*abs(q*cd-sd)/(q*W);
            X=Xnext;
        end
        scale=(o.t(j+1)/t)^2;D=diag([scale,1]);X=D*X*D;
        o.A(j+1)=X(1,1);
    end
    o.r=exp([flipud(cumsum(flipud(o.logW)));0]);
    o.theta=[0;cumsum(o.delta)];
    o.x=o.r.*[cos(o.theta),sin(o.theta)];o.g=o.zeta.*o.x;
    tail=max(1,N-9999):N;
    o.summary=struct('J',J,'cycles',N,'nu_over_a',nuRatio,'initial_scale',t0, ...
        'turns',sum(o.delta)/(2*pi),'initial_radius',o.r(1),'terminal_radius',o.r(end), ...
        'initial_gradient_norm',norm(o.g(1,:)),'final_gradient_norm',norm(o.g(end,:)), ...
        'median_A',median(o.A(tail)),'target_A',K, ...
        'median_QE',median(o.Q(2*tail-1)),'target_QE',z+2*a, ...
        'median_QO',median(o.Q(2*tail)),'target_QO',w+a, ...
        'min_armijo',min(o.armijo),'max_strong',max(o.strong), ...
        'median_strong_even',median(o.strong(2*tail-1)), ...
        'median_strong_odd',median(o.strong(2*tail)), ...
        'target_strong_even',1-2*a/(z+2*a),'target_strong_odd',1-a/(w+a), ...
        'max_secant_residual',max(o.secant),'max_update_residual',max(o.update), ...
        'armijo_failures',sum(o.armijo<.25-1e-12),'strong_failures',sum(o.strong>.75+1e-12), ...
        'min_normalized_eigenvalue',o.minEig,'min_sigma',o.minSigma);
end
function obj=finite(o)
    obj.points=o.x;obj.corrections=(o.zeta-1).*o.x;
    obj.tree=KDTreeSearcher(obj.points);
    [~,d]=knnsearch(obj.tree,obj.points,'K',2);
    obj.radii=.24*d(:,2);assert(all(obj.radii>0));
    obj.ratio=max(vecnorm(obj.corrections,2,2)./obj.radii);
    obj.perturbation=(135/16+15*sqrt(3)/2)*obj.ratio;
    obj.band=[1-obj.perturbation,1+obj.perturbation];
    obj.steps=2*o.cycles;obj.J=o.J;
end
function [f,g]=valueGrad(obj,x)
    f=.5*(x'*x);g=x;
    [i,d]=knnsearch(obj.tree,x');rho=obj.radii(i);
    if d>=rho,return,end
    v=x-obj.points(i,:)';corr=obj.corrections(i,:)';u=d/rho;
    if u<=1/3,phi=1;dp=0;else
        z=(3*u-1)/2;phi=1-10*z^3+15*z^4-6*z^5;
        dp=-45*z^2*(1-z)^2;
    end
    cv=corr'*v;f=f+phi*cv;g=g+phi*corr;
    if d>0 && dp~=0,g=g+cv*dp*v/(rho*d);end
end
function result=solve(obj,o,method,maxIter)
    x0=o.x(1,:)';L=chol(o.h0,'lower');tol=1e-10;
    trace=zeros(maxIter+1,5);nr=0;neval=0;lastZ=[0;0];caught='';
    opts=optimset('Algorithm','quasi-newton','HessUpdate',method,'GradObj','on', ...
        'Display','off','TolFun',0,'TolX',0,'MaxIter',maxIter,'MaxFunEvals',100000,'OutputFcn',@observe);
    started=tic;
    try
        [zz,fval,flag,out]=fminunc(@evaluate,[0;0],opts);
    catch err
        if ~strcmp(err.identifier,'optim:lineSearch:FPrimeInitialNeg'),rethrow(err);end
        zz=lastZ;[fval,~]=evaluate(zz);flag=NaN;caught=err.message;
        out=struct('iterations',trace(max(1,nr),1),'funcCount',neval,'message',err.message,'identifier',err.identifier);
    end
    elapsed=toc(started);xx=x0+L*zz;[fv,gg]=PolarDFP.valueGrad(obj,xx);
    assert(abs(fv-fval)<1e-12*max(1,abs(fval)));
    if ~isempty(caught),status='native_solver_error';elseif norm(gg)<=tol,status='gradient_tolerance';elseif out.iterations>=maxIter,status='iteration_limit';else,status='solver_stopped';end
    result.trace=array2table(trace(1:nr,:),'VariableNames',{'iteration','function_value','gradient_norm','x1','x2'});
    result.x=xx;result.output=out;
    result.summary=struct('method',method,'J',o.J,'prefix_steps',2*o.cycles, ...
        'status',status,'exitflag',flag,'iterations',out.iterations,'function_evaluations',out.funcCount, ...
        'final_gradient_norm',norm(gg),'elapsed_seconds',elapsed,'hessian_band',obj.band, ...
        'max_iterations',maxIter,'gradient_tolerance',tol,'matlab_version',version,'solver_message',out.message);
    function [f,gz]=evaluate(zz)
        [f,gx]=PolarDFP.valueGrad(obj,x0+L*zz);gz=L'*gx;neval=neval+1;
    end
    function stop=observe(zz,values,state)
        lastZ=zz;xx=x0+L*zz;[ff,gg]=PolarDFP.valueGrad(obj,xx);
        stop=norm(gg)<=tol||values.iteration>=maxIter;
        if strcmp(state,'iter')||strcmp(state,'init')
            nr=nr+1;trace(nr,:)=[values.iteration,ff,norm(gg),xx'];
        end
    end
end
function writeJSON(path,s)
    fid=fopen(path,'w');assert(fid>0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(s,'PrettyPrint',true));
end
end
end
