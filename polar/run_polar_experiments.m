function run_polar_experiments(mode, outputDir)
% SPDX-License-Identifier: MIT
% RUN_POLAR_EXPERIMENTS Reproduce the prescribed polar DFP experiments.
%   run_polar_experiments('all') runs both comparisons and exact certificates.
%   'scales' runs only the three beta recurrence diagnostics and their figure.
%   'smoke', 'wolfe', 'certify', and 'figures' select the corresponding stages.
if nargin<1,mode='all';end
mode=validatestring(mode,{'all','smoke','wolfe','scales','certify','figures'});
sourceDir=fileparts(mfilename('fullpath'));
if nargin<2,outputDir=fullfile(sourceDir,'results');end
if ~isfolder(outputDir),mkdir(outputDir);end
previousPath=path;restorePath=onCleanup(@()path(previousPath));
addpath(sourceDir);previous=pwd;cd(outputDir);cleanup=onCleanup(@()cd(previous));
if strcmp(mode,'all')
    for stale={'complete.json','certificate.json','finite_200_certificate.json'}
        if isfile(stale{1}),delete(stale{1});end
    end
end
if strcmp(mode,'smoke')&&isfile('smoke.json'),delete('smoke.json');end
if any(strcmp(mode,{'all','figures'}))&&isfile('figures_complete.json'),delete('figures_complete.json');end
if strcmp(mode,'figures'),make_polar_figures;return,end
if strcmp(mode,'certify'),certify_finite;certify_finite('finite_200');return,end
if any(strcmp(mode,{'all','wolfe'}))
    pairs=[.25,.75;.2,.3;.8,.9];
    for caseIndex=1:size(pairs,1)
        if caseIndex==1
            base=2*pi/9;
            pars=PolarDFP.parameters(2*base,base,1.03*base,.5,pairs(caseIndex,1),pairs(caseIndex,2));
            stem='wolfe_reference';
        else
            pars=PolarDFP.wolfeParameters(pairs(caseIndex,1),pairs(caseIndex,2),.5);
            stem=sprintf('wolfe_case_%d',caseIndex-1);
        end
        diagnostic=PolarDFP.orbit(10000,400^3,pars);
        assert(diagnostic.summary.armijo_failures==0 && diagnostic.summary.strong_failures==0);
        assert(diagnostic.summary.max_update_residual<1e-10);
        PolarDFP.writeJSON([stem,'.json'],diagnostic.summary);
        save([stem,'.mat'],'diagnostic','-v7.3');
    end
    if strcmp(mode,'wolfe'),return,end
end
if any(strcmp(mode,{'all','scales'}))
    % Vary beta alone, retaining the representative geometric shape and Wolfe pair.
    betas=[1/3,1/2,2/3];base=2*pi/9;
    for caseIndex=1:numel(betas)
        pars=PolarDFP.parameters(2*base,base,1.03*base,betas(caseIndex),.25,.75);
        diagnostic=scaling_diagnostics(10000000,32^3,pars);
        assert(diagnostic.summary.armijo_failures==0 && diagnostic.summary.strong_failures==0);
        assert(diagnostic.summary.sampled_max_update_residual<1e-10);
        assert(all(isfinite(diagnostic.ratio_radial)) && all(diagnostic.ratio_radial>0));
        assert(all(isfinite(diagnostic.ratio_correction)) && all(diagnostic.ratio_correction>0));
        stem=sprintf('exponent_case_%d',caseIndex);
        PolarDFP.writeJSON([stem,'.json'],diagnostic.summary);
        save([stem,'.mat'],'diagnostic','-v7.3');
    end
    PolarDFP.writeJSON('scales_environment.json',struct('matlab',version, ...
        'architecture',computer,'threads',maxNumCompThreads,'precision','IEEE 754 binary64'));
    if strcmp(mode,'scales'),make_scaling_figure;return,end
end
if strcmp(mode,'smoke')
    o=PolarDFP.orbit(200,400^3,1.03);obj=PolarDFP.finite(o);
    r=PolarDFP.solve(obj,o,'bfgs',100);
    assert(o.summary.armijo_failures==0&&o.summary.strong_failures==0);
    assert(obj.band(1)>0&&r.summary.final_gradient_norm<=1e-10);
    PolarDFP.writeJSON('smoke.json',struct('recurrence',o.summary,'solver',r.summary));return
end
env=struct('matlab',version,'architecture',computer,'threads',maxNumCompThreads, ...
    'precision','IEEE 754 binary64','optimization_toolbox',ver('optim'), ...
    'statistics_toolbox',ver('stats'),'symbolic_toolbox',ver('symbolic'));
if ismac,[~,cpu]=system('sysctl -n machdep.cpu.brand_string');env.processor=strtrim(cpu);end
PolarDFP.writeJSON('environment.json',env);
geometry=PolarDFP.orbit(100000,32^3,1.03);
PolarDFP.writeJSON('geometry.json',geometry.summary);save('geometry.mat','geometry','-v7.3');
for rootJ=[400,200]
    if rootJ==400,suffix='';else,suffix='_200';end
    finiteOrbit=PolarDFP.orbit(4002,rootJ^3,1.03);objective=PolarDFP.finite(finiteOrbit);
    stem=['finite',suffix];save([stem,'.mat'],'finiteOrbit','objective','-v7.3');
    writematrix([objective.points,objective.corrections,objective.radii],[stem,'_fixed_data.csv']);
    PolarDFP.writeJSON([stem,'.json'],struct('orbit',finiteOrbit.summary, ...
        'hessian_band',objective.band,'ratio',objective.ratio,'points',size(objective.points,1)));
    dfp=PolarDFP.solve(objective,finiteOrbit,'dfp',5000);
    bfgs=PolarDFP.solve(objective,finiteOrbit,'bfgs',5000);
    PolarDFP.writeJSON(['dfp',suffix,'.json'],dfp.summary);
    PolarDFP.writeJSON(['bfgs',suffix,'.json'],bfgs.summary);
    writetable(dfp.trace,['dfp',suffix,'.csv']);writetable(bfgs.trace,['bfgs',suffix,'.csv']);
    save(['dfp',suffix,'.mat'],'dfp');save(['bfgs',suffix,'.mat'],'bfgs');
    certify_finite(stem);
end
make_polar_figures;
PolarDFP.writeJSON('complete.json',struct('status','completed','matlab',version));
end
