classdef PolarDFP
% SPDX-License-Identifier: MIT
% Exact prescribed polar angles and finite interpolation for classical DFP.
methods (Static)
function pars=parameters(a,b,nu,beta,c1,c2)
    % General two-stage cycle; a and b are the even/odd angle coefficients.
    assert(beta>0 && beta<1 && 0<c1 && c1<c2 && c2<1);
    assert(a>0 && b>0 && nu>max(0,(a-b)/2) && nu<a);
    den=b-a+2*nu;
    le=nu*(b+nu)/den;lo=nu*(a-nu)/den;
    gamma=den/(nu*(a-nu)*(b+nu));
    mult=(1-1/(gamma*le^2))*(1-1/(gamma*lo^2));
    we=lo/(a+lo);wo=le/(b+le);
    assert(abs(mult)<1,'Unstable transverse two-cycle');
    assert(max([we,wo])<c2 && min([(1+we)/2,(1+wo)/2])>c1, ...
        'No strict limiting Wolfe margins');
    pars=struct('a',a,'b',b,'nu',nu,'beta',beta,'p',1/(1+beta), ...
        'c1',c1,'c2',c2,'lambdaE',le,'lambdaO',lo,'gamma',gamma, ...
        'transverse_multiplier',mult,'strong_even',we,'strong_odd',wo);
end
function pars=wolfeParameters(c1,c2,beta)
    % Construct one stable shape for any fixed standard Wolfe pair.
    if nargin<3,beta=.5;end
    assert(0<c1 && c1<c2 && c2<1 && beta>0 && beta<1);
    if c2<=1/3
        ell=max(4,1/c2+.5);z=ell;
    else
        w=(max(1/3,2*c1-1)+c2)/2;z=1/w-1;
        found=false;
        for k=1:52
            ell=1+2^(-k);a=ell-1+z;b=1-ell+z*ell;
            we=1/(ell+z);wo=ell/(1+z*ell);m=(1-z/ell)*(1-z*ell);
            if b>0 && abs(m)<1 && max(we,wo)<c2 && min((1+we)/2,(1+wo)/2)>c1
                found=true;break
            end
        end
        assert(found,'Wolfe margins are below binary64 parameter resolution');
    end
    pars=PolarDFP.parameters(ell-1+z,1-ell+z*ell,ell-1,beta,c1,c2);
    pars.ell=ell;pars.z=z;
end
function nu=inverseGamma(gamma,pars,nu)
    % Follow the local increasing inverse branch used in the proof.
    a=pars.a;b=pars.b;lower=max(0,(a-b)/2);
    for k=1:12
        value=(b-a+2*nu)/(nu*(a-nu)*(b+nu));
        if abs(value-gamma)<=8*eps*max(abs(gamma),abs(value)),return,end
        deriv=value*(2/(b-a+2*nu)-1/nu+1/(a-nu)-1/(b+nu));
        assert(deriv>0 && isfinite(deriv),'Left the local increasing inverse branch');
        change=(value-gamma)/deriv;
        change=max(-.5*(a-nu),min(.5*(nu-lower),change));
        nu=nu-change;
    end
    value=(b-a+2*nu)/(nu*(a-nu)*(b+nu));
    assert(abs(value-gamma)<=1e-12*abs(gamma),'Inverse-branch solve did not converge');
end
function o=orbit(N,J,nuRatio)
    % A numeric third argument retains the paper's 2:1 representative shape.
    % nuRatio denotes nu_star/b, including in the saved metadata.
    % Pass PolarDFP.wolfeParameters(c1,c2,beta) for a general Wolfe pair.
    if nargin<3,nuRatio=1.03;end
    validateattributes(N,{'numeric'},{'scalar','integer','positive'});
    validateattributes(J,{'numeric'},{'scalar','positive','finite'});
    if isstruct(nuRatio)
        pars=nuRatio;nuRatio=pars.nu/pars.b;
    else
        base=2*pi/9;
        pars=PolarDFP.parameters(2*base,base,nuRatio*base,.5,.25,.75);
    end
    a=pars.a;b=pars.b;p=pars.p;nu=pars.nu;
    z=pars.lambdaO;w=pars.lambdaE;K=pars.gamma;
    phi0=p/(nu*(1/a+1/b));X=[K,K*w;K*w,1+K*w*w];
    n=2*N;tau0=J^(-p);S=diag([tau0,1]);
    o.h0=S*(X\S);o.h0=(o.h0+o.h0')/2;
    o.J=J;o.nuRatio=nuRatio;o.cycles=N;o.parameters=pars;
    o.feedback='phi=p*gamma_star/(gamma*nu(gamma)*(1/a+1/b))';
    o.targets=struct('A',K,'QE',z+a,'QO',w+b,'D0',a*(z+a)-a*a/2,'D1',b*(w+b)-b*b/2,'L',phi0);
    o.tau=(J+(0:N)') .^(-p);o.t=sqrt(o.tau);o.A=zeros(N+1,1);o.A(1)=K;
    o.phi=zeros(N,1);o.nu=zeros(N,1);
    o.logW=zeros(n,1);o.delta=zeros(n,1);o.Q=zeros(n,1);
    o.armijo=zeros(n,1);o.strong=zeros(n,1);o.secant=zeros(n,1);o.update=zeros(n,1);
    o.zeta=ones(n+1,1);o.minEig=Inf;o.minSigma=Inf;
    for j=1:N
        tau=o.tau(j);index=J+j-1;v=1/index;
        angleFactor=expm1((1-p)*log1p(v))/((1-p)*v);
        nu=PolarDFP.inverseGamma(X(1,1),pars,nu);
        phi=p*K/(X(1,1)*nu*(1/a+1/b));o.phi(j)=phi;o.nu(j)=nu;
        dz=expm1(-phi*v);
        factors=[1,1+dz,1];o.zeta(2*j)=1+dz;
        for stage=1:2
            k=2*j-2+stage;coeff=[a,b];delta=coeff(stage)*tau*angleFactor;
            zet=factors(stage);zhat=factors(stage+1);
            if stage==1,change=dz;else,change=-dz;end
            Q=X(2,2)/X(1,2);q=tau*Q;
            vb=[-Q;1];wb=[tau*change*cot(delta)-zet*tau*tau*Q;zhat];
            sigma=vb'*wb;P=eye(2)-wb*vb'/sigma;
            Y=P*X*P'+wb*wb'/sigma;Y=(Y+Y')/2;
            sd=sin(delta);cd=cos(delta);Rbar=[cd,tau*sd;-sd/tau,cd];
            Xnext=Rbar*Y*Rbar';Xnext=(Xnext+Xnext')/2;
            o.minEig=min(o.minEig,min(eig(Xnext)));o.minSigma=min(o.minSigma,sigma);
            assert(sigma>0 && Xnext(1,2)>0 && det(Xnext)>0,'Invalid polar state');
            vv=[-q;1];ww=[wb(1)/tau;wb(2)];
            H=diag([tau,1])*(X\diag([tau,1]));
            Hnew=diag([tau,1])*(Y\diag([tau,1]));Hy=H*ww;
            Hdfp=H-Hy*Hy'/(ww'*Hy)+vv*vv'/sigma;
            o.secant(k)=norm(Hnew*ww-vv)/norm(vv);
            o.update(k)=norm(Hdfp-Hnew,'fro')/norm(Hnew,'fro');
            wm1=-2*sin(delta/2)^2+q*sd;
            assert(wm1>0,'Radius did not decrease');W=1+wm1;
            o.logW(k)=log1p(wm1);o.delta(k)=delta;o.Q(k)=Q;
            o.armijo(k)=wm1*(2+wm1)/(2*zet*q*sd*W);
            o.strong(k)=zhat/zet*abs(q*cd-sd)/(q*W);X=Xnext;
        end
        scale=o.tau(j+1)/tau;D=diag([scale,1]);X=D*X*D;o.A(j+1)=X(1,1);
    end
    o.r=exp([flipud(cumsum(flipud(o.logW)));0]);o.theta=[0;cumsum(o.delta)];
    o.x=o.r.*[cos(o.theta),sin(o.theta)];o.g=o.zeta.*o.x;
    % Dimensionless diagnostics use stable differences from the recurrence.
    % The cycle index is n_j=J+j, with j=0,...,N-1.
    cycleIndex=J+(0:N-1)';
    o.ratio_radial=(o.logW(1:2:end)+o.logW(2:2:end))./ ...
        ((o.targets.D0+o.targets.D1)*o.tau(1:end-1).^2);
    o.ratio_correction=cycleIndex.*(-expm1(-o.phi./cycleIndex))/o.targets.L;
    tail=max(1,N-9999):N;
    o.summary=struct('J',J,'cycles',N,'nu_over_b',nuRatio,'initial_scale',sqrt(tau0), ...
        'beta',pars.beta,'p',p,'angle_even',a,'angle_odd',b,'nu_star',pars.nu, ...
        'wolfe_c1',pars.c1,'wolfe_c2',pars.c2,'feedback',o.feedback, ...
        'turns',sum(o.delta)/(2*pi),'initial_radius',o.r(1),'terminal_radius',o.r(end), ...
        'initial_gradient_norm',norm(o.g(1,:)),'final_gradient_norm',norm(o.g(end,:)), ...
        'median_A',median(o.A(tail)),'target_A',K, ...
        'median_QE',median(o.Q(2*tail-1)),'target_QE',z+a, ...
        'median_QO',median(o.Q(2*tail)),'target_QO',w+b, ...
        'target_D0',o.targets.D0,'target_D1',o.targets.D1,'target_L',o.targets.L, ...
        'min_ratio_radial',min(o.ratio_radial),'max_ratio_radial',max(o.ratio_radial), ...
        'median_ratio_radial',median(o.ratio_radial(tail)),'final_ratio_radial',o.ratio_radial(end), ...
        'min_ratio_correction',min(o.ratio_correction),'max_ratio_correction',max(o.ratio_correction), ...
        'median_ratio_correction',median(o.ratio_correction(tail)), ...
        'final_ratio_correction',o.ratio_correction(end), ...
        'min_armijo',min(o.armijo),'max_strong',max(o.strong), ...
        'median_strong_even',median(o.strong(2*tail-1)), ...
        'median_strong_odd',median(o.strong(2*tail)), ...
        'target_strong_even',pars.strong_even,'target_strong_odd',pars.strong_odd, ...
        'max_secant_residual',max(o.secant),'max_update_residual',max(o.update), ...
        'armijo_failures',sum(o.armijo<pars.c1-1e-12),'strong_failures',sum(o.strong>pars.c2+1e-12), ...
        'min_normalized_eigenvalue',o.minEig,'min_sigma',o.minSigma, ...
        'min_phi',min(o.phi),'max_phi',max(o.phi),'transverse_multiplier',pars.transverse_multiplier);
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
        zz=lastZ;fval=PolarDFP.valueGrad(obj,x0+L*zz);flag=NaN;caught=err.message;
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
            if nr==0||trace(nr,1)~=values.iteration,nr=nr+1;end
            trace(nr,:)=[values.iteration,ff,norm(gg),xx'];
        end
    end
end
function writeJSON(path,s)
    fid=fopen(path,'w');assert(fid>0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(s,'PrettyPrint',true));
end
end
end
