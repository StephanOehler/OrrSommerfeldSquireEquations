% This file generates the support files for Figure 11 of
% Int J Heat Fluid Flow (2021), vol. 87, pp. 108735

clear


alphas = [1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2];
names = {'n5','n4','n3','n2','n1','1','p1','p2'};


IOC_vec_a = size(alphas);
ME_vec = size(alphas);



for i = 1:length(alphas)
    
    load(['..\Week_260_Reviewers_FIX\Cost_change/alpha_',names{i}],'gam_z','gam_uz')

    
    IOC_vec_a(i) = gam_uz;
    ME_vec(i) = gam_z;
    
end

compounded_results1 = real([ME_vec;IOC_vec_a]).';

figure('position',[0,0,2*sqrt(2),1]*200,'color','white')
bb = bar([1:length(alphas)],compounded_results1);


set(gca,'TickLabelInterpreter','Latex')
set(gca,'XTickLabel',{'$10^{-5}$','$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$1$','$10^{1}$','$10^{2}$'})

ylabel('$\mathbf{E}_{u,v,w}$','interpreter','latex')
xlabel('$\alpha$','interpreter','latex')


bb(1).FaceColor = [1,1,1] * .7;
bb(2).FaceColor = [1,1,1] * 0.3;


print('fig11','-depsc')