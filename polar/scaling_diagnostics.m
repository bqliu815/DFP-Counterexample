function d = scaling_diagnostics(N,J,pars)
% SCALING_DIAGNOSTICS Run the polar recurrence with bounded output storage.
% Every cycle is evaluated. Only plotting samples and expensive residual
% checks are thinned; Wolfe ratios and matrix validity are checked at every step.
% SPDX-License-Identifier: MIT
validateattributes(N,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(J,{'numeric'},{'scalar','positive','finite'});
a=pars.a;b=pars.b;p=pars.p;nu=pars.nu;K=pars.gamma;w=pars.lambdaE;z=pars.lambdaO;
phi0=p/(nu*(1/a+1/b));
D0=a*(z+a)-a*a/2;D1=b*(w+b)-b*b/2;Dsum=D0+D1;
X=[K,K*w;K*w,1+K*w*w];
cycles=unique([0:min(100,N-1),round(logspace(2,log10(max(100,N-1)),10000)), ...
    99999,999999,9999999,N-1])';
cycles=cycles(cycles<N);count=numel(cycles);nextSample=1;
d=struct('J',J,'cycles',N,'parameters',pars,'cycle',cycles, ...
    'ratio_radial',zeros(count,1),'ratio_correction',zeros(count,1), ...
    'gamma',zeros(count,1),'phi',zeros(count,1));
minSigma=Inf;minDet=Inf;minEig=Inf;minArmijo=Inf;maxStrong=0;
maxSecant=0;maxUpdate=0;armijoFailures=0;strongFailures=0;
maxRadial=-Inf;maxCorrection=-Inf;peakRadial=0;peakCorrection=0;
started=tic;
% Small blocks retain the same vectorized powers as PolarDFP.orbit.
for first=1:100000:N
    last=min(first+99999,N);tauBlock=(J+(first-1:last)').^(-p);
    for j=first:last
        pos=j-first+1;tau=tauBlock(pos);index=J+j-1;v=1/index;
        angleFactor=expm1((1-p)*log1p(v))/((1-p)*v);
        nu=PolarDFP.inverseGamma(X(1,1),pars,nu);
        phi=p*K/(X(1,1)*nu*(1/a+1/b));dz=expm1(-phi*v);
        factors=[1,1+dz,1];logW=zeros(1,2);gamma=X(1,1);
        sampleNow=nextSample<=count && j-1==cycles(nextSample);
        for stage=1:2
            coeff=[a,b];delta=coeff(stage)*tau*angleFactor;
            zet=factors(stage);zhat=factors(stage+1);
            if stage==1,change=dz;else,change=-dz;end
            Q=X(2,2)/X(1,2);q=tau*Q;
            vb=[-Q;1];wb=[tau*change*cot(delta)-zet*tau*tau*Q;zhat];
            sigma=vb'*wb;P=eye(2)-wb*vb'/sigma;
            Y=P*X*P'+wb*wb'/sigma;Y=(Y+Y')/2;
            sd=sin(delta);cd=cos(delta);Rbar=[cd,tau*sd;-sd/tau,cd];
            Xnext=Rbar*Y*Rbar';Xnext=(Xnext+Xnext')/2;
            determinant=det(Xnext);
            assert(sigma>0 && Xnext(1,1)>0 && Xnext(1,2)>0 && determinant>0 ...
                && all(isfinite(Xnext(:))),'Invalid polar state');
            minSigma=min(minSigma,sigma);minDet=min(minDet,determinant);
            if sampleNow
                minEig=min(minEig,min(eig(Xnext)));
                vv=[-q;1];ww=[wb(1)/tau;wb(2)];
                H=diag([tau,1])*(X\diag([tau,1]));
                Hnew=diag([tau,1])*(Y\diag([tau,1]));Hy=H*ww;
                Hdfp=H-Hy*Hy'/(ww'*Hy)+vv*vv'/sigma;
                maxSecant=max(maxSecant,norm(Hnew*ww-vv)/norm(vv));
                maxUpdate=max(maxUpdate,norm(Hdfp-Hnew,'fro')/norm(Hnew,'fro'));
            end
            wm1=-2*sin(delta/2)^2+q*sd;
            assert(wm1>0 && isfinite(wm1),'Radius did not decrease');W=1+wm1;
            logW(stage)=log1p(wm1);
            armijo=wm1*(2+wm1)/(2*zet*q*sd*W);
            strong=zhat/zet*abs(q*cd-sd)/(q*W);
            minArmijo=min(minArmijo,armijo);maxStrong=max(maxStrong,strong);
            armijoFailures=armijoFailures+(armijo<pars.c1-1e-12);
            strongFailures=strongFailures+(strong>pars.c2+1e-12);X=Xnext;
        end
        scale=tauBlock(pos+1)/tau;D=diag([scale,1]);X=D*X*D;
        radial=(logW(1)+logW(2))/(Dsum*tau^2);
        correction=index*(-expm1(-phi/index))/phi0;
        assert(all(isfinite([radial,correction])));
        if radial>maxRadial,maxRadial=radial;peakRadial=j-1;end
        if correction>maxCorrection,maxCorrection=correction;peakCorrection=j-1;end
        if sampleNow
            d.ratio_radial(nextSample)=radial;d.ratio_correction(nextSample)=correction;
            d.gamma(nextSample)=gamma;d.phi(nextSample)=phi;nextSample=nextSample+1;
        end
    end
end
d.summary=struct('J',J,'cycles',N,'beta',pars.beta,'p',p, ...
    'angle_even',a,'angle_odd',b,'nu_star',pars.nu,'wolfe_c1',pars.c1,'wolfe_c2',pars.c2, ...
    'stored_cycles',count,'target_radial',1,'target_correction',1, ...
    'max_ratio_radial',maxRadial,'peak_radial_cycle',peakRadial, ...
    'max_ratio_correction',maxCorrection,'peak_correction_cycle',peakCorrection, ...
    'final_ratio_radial',radial,'final_ratio_correction',correction,'final_gamma',X(1,1), ...
    'min_armijo',minArmijo,'max_strong',maxStrong, ...
    'armijo_failures',armijoFailures,'strong_failures',strongFailures, ...
    'min_sigma',minSigma,'min_normalized_determinant',minDet, ...
    'sampled_min_normalized_eigenvalue',minEig, ...
    'sampled_max_secant_residual',maxSecant,'sampled_max_update_residual',maxUpdate, ...
    'elapsed_seconds',toc(started));
end
