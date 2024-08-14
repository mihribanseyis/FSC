%22/24 single gene knock-outs:
clear
KOIDs = ["b0756" "b2388" "b0688" "b4025" "b3916" "b1723" "b4232" "b2097" "b0755" "b1854" "b1676" "b1702" "b1852" "b0767" "b2029" "b3386" "b2914" "b4090" "b2935" "b2465" "b2464" "b0008"];

load("iML1515.mat");
model = iML1515;
pattern = 'b\d{4}';

%ub/lb depending on growth rates
%     for   0.2   0.1  0.4  0.5  0.7
growth_ubs=[0.22 0.12 0.42 0.52 0.72];
growth_lbs=[0.18 0.08 0.38 0.48 0.68];

for k = 1:numel(KOIDs)
    KOgene = KOIDs(k); % current gene KO
    for i = 1:size(model.grRules)
        grRule = model.grRules{i}; % grRule for current reaction
        model.KOs.lb(k,i) = model.lb(i);
        model.KOs.ub(k,i) = model.ub(i);
        % Check if the current gene is present in the current reaction's rule
        if contains(grRule, KOgene)
            orConditions = split(grRule, ' or ');
            grCheck = zeros(1, numel(orConditions));
            for oridx=1:size(orConditions)
                andConditions = split(orConditions{oridx}, ' and ');
                andCheck = zeros(1, numel(andConditions));
                for andidx=1:size(andConditions)
                    matches = regexp(andConditions{andidx}, pattern, 'match'); %Extract gene names
                    KOcheck = contains(KOgene, matches{1}); %KO or not
                    if ~KOcheck
                        andCheck(andidx) = 1; % GPR rule has no KO
                    else
                        andCheck(andidx) = 0; % GPR rule has KO
                        break; % Exit the loop if any KO is found within an 'and' group
                    end
                end
                if sum(andCheck) ~= numel(andConditions) %all KOs in 'and' condition
                    grCheck(oridx) = 0;
                else
                    grCheck(oridx) = 1;
                end
            end
            % Apply knockout only for this condition
            if sum(grCheck) == 0 %all KOs in or cond.
                model.KOs.lb(k, i) = 0; % Lower bound set to 0 (knockout)
                model.KOs.ub(k, i) = 0; % Upper bound set to 0 (knockout)
            end
        end
    end
end

%pFBA application - LP
%opt biomass & weight according to KOs
biomass=find(model.c==1);
[m,r]=size(model.S);
weigth=zeros(22,r);

for wix=1:22
    for wiy=1:r
        if model.KOs.ub(wix, wiy) ~= 0
            weigth(wix,wiy) = 1;
        end
    end
end

%inequality constraints
A=[];
b=[];

%equality constraints
%first eq -> N.V=0
%second eq -> V-(Vi+)+(Vi-)=0

Aeq_1=[model.S zeros(m,r) zeros(m,r)];
Aeq_2=[eye(r) -eye(r) eye(r)];
Aeq_3=zeros(1,3*r);
Aeq=[Aeq_1;Aeq_2];

beq_1=zeros(m,1);
beq_2=zeros(r,1);
beq=[beq_1;beq_2];

%objective function
f=[zeros(1,r) weigth(1,:) weigth(1,:)];

%flux distributions for different growth rates
for gr=1:5
    model.ub(biomass)=growth_ubs(gr);
    model.lb(biomass)=growth_lbs(gr);
    ub=[model.ub;model.ub;model.ub];
    lb=[model.lb;zeros(r,1);zeros(r,1)];
    [sol.x(gr,:),sol.fval(gr,:)]=linprog(f,A,b,Aeq,beq,lb,ub);
end

%flux distributions for KOs
model.ub(biomass)=growth_ubs(1);
model.lb(biomass)=growth_lbs(1);
for ko_flx=1:22
    ub=[model.KOs.ub(ko_flx,:).';model.KOs.ub(ko_flx,:).';model.KOs.ub(ko_flx,:).'];
    lb=[model.KOs.lb(ko_flx,:).';zeros(r,1);zeros(r,1)];
    f=[zeros(1,r) weigth(ko_flx,:) weigth(ko_flx,:)];
    [sol.x(5+ko_flx,:),sol.fval(5+ko_flx,:)]=linprog(f,A,b,Aeq,beq,lb,ub);
end

%flux sums
%phi=0.5*sum(|Sij*vj|)
for condi=1:27
    for meti=1:m
        sol.flux_sums(condi,meti) = FluxSum(model.S(meti,:), sol.x(condi,1:r));
    end
end

sol.flux_sums=sol.flux_sums';

%metabolites from Ishii
%find(strcmp(model.mets, 'dhap_c'))
met_EIs = readcell("met_EIs.xlsx");
%col_EI = size(met_EIs,2);
for isi=3:size(met_EIs)
    ishii.mets(isi-2,1) = met_EIs(isi,3);
end

%EI values of 114 mets
corrl.EIs = cell(114, 27);
for eix=1:114
    for eiy=1:27
        corrl.EIs(eix,eiy) = met_EIs(eix+2,eiy+9);
    end
end

%extracting the flux sums for measured metabolites
%removing the _c,_e,_p from met names
corrl.mets = cell(m, 1);
for fli=1:m
    met_cur = split(model.mets(fli),'_');
    corrl.mets{fli} = met_cur{1};
end

corrl.sums=zeros(114,27);
corrl.sums=[];
for im=1:114
    %indexes of mets in model that are both in ishii&model
    met_idx = find(strcmp(corrl.mets, ishii.mets(im)));
    if ~isempty(met_idx)
        if size(met_idx,1) == 1
            for mti=1:27
                corrl.sums(im,mti)=sol.flux_sums(met_idx(1), mti);
            end
        else
            mtsms=zeros(1,27);
            for ix=1:numel(met_idx)
                for mt2i=1:27
                    mtsms(mt2i)=mtsms(mt2i)+sol.flux_sums(met_idx(ix),mt2i);
                    corrl.sums(im,mt2i) = mtsms(mt2i);
                end
            end
        end
    end
end

%removing the 'NA' values
EI_NA=readtable("met_EIs.xlsx",Sheet="Sheet2");
EI_NA=EI_NA(:,1:end-1);
missing_indices=ismissing(EI_NA);
%corrl.sums(missing_indices)=NaN;

%how many conditions for each metabolite is reported
corrl.NAs = zeros(114,1);
for ina=1:114
    corrl.NAs(ina) = 27-sum(ismissing(EI_NA(ina,:)));
end

%correlation
EI_NA= table2array(EI_NA);
corrl.res.rho = zeros(114,1);
corrl.res.pval = zeros(114,1);

for row=1:114
    row_1 = EI_NA(row,:);
    row_2 = corrl.sums(row,:);
    %non-NaN indices
    nonNan_row1 = ~isnan(row_1);
    %nonNan_row2 = ~isnan(row_2);
    row1_filtered = row_1(nonNan_row1);
    row2_filtered = row_2(nonNan_row1);
    [corrl.res.rho(row,1),corrl.res.pval(row,1)] = corr(row1_filtered', row2_filtered');
end

filename = 'met_EIs.xlsx';

function phi = FluxSum(S, vj)
    phi = 0;
    for idxi=1:2712
        phi=phi+abs(S(idxi)*vj(idxi));
    end
    phi = phi*0.5;
end