%initCobraToolbox(false)
%FScoupling function
function [mins maxs couplings] = FScoupling(model)
%model=convertToIrreversible(model);
[m,r]=size(model.S);

mins = zeros(m,m);
maxs = zeros(m,m);
couplings = zeros(m,m);

%var. vector -> [v1' v2' ... vr' t]'
%Equality constraints
%Sv'=0 & d^t.v_j'=1
Aeq = [model.S zeros(m,1); zeros(1,r+1)];
beq = [zeros(m+1,1)];
beq(end)=1;

%inequality constraints
%v'-ub.t <= 0 & lb.t-v' <= 0
A = [eye(r) -model.ub; -eye(r) model.lb];
b = [zeros(2*r,1)];

%ub&lb
ub = [model.ub; 1000];
lb = [model.lb; -1000];
%80% biomass
%biomass=find(model.c==1);
%lb(biomass)=model.ub(biomass)*0.1;

for i=1:m
    for j=1:m
        if i~=j
            %inequality constraints are the same for each couple
            %equality constraints - d^T for each j is different
            Aeq(end,1:r) = abs(model.S(j,:));
            %objective function
            f = zeros(r+1,1);
            %minimization 
            f(1:r)= abs(model.S(i,:)); 
            [x fval] = linprog(f,A,b,Aeq,beq,lb,ub);
            if ~isempty(fval)
                mins(i,j)= fval;
            else 
                mins(i,j) = NaN;
            end
            %maximization
            f(1:r)= -abs(model.S(i,:));
            [x fval] = linprog(f,A,b,Aeq,beq,lb,ub);
            if ~isempty(fval)
                maxs(i,j) = -fval;
            else
                maxs(i,j) = NaN;
            end
        end
    end
end

%coupling classification
for ix=1:m
    for jx=1:m
        if ix~=jx
            % uncoupled - 0
            if mins(ix,jx) == 0 && maxs(ix,jx)>10000
                couplings(ix,jx) = 0;
            % fully coupled - 1
            %elseif mins(ix,jx) == maxs(ix,jx) && maxs(ix,jx)>0
            elseif abs(mins(ix,jx)-maxs(ix,jx))<0.001 && mins(ix,jx)>0
                couplings(ix,jx) = 1;
            % partially coupled - 2  maxs ~= +Inf ?
            elseif mins(ix,jx)>0 && maxs(ix,jx)>0 && mins(ix,jx)~=maxs(ix,jx)
                couplings(ix,jx) = 2;
            % directionally (i->j) coupled - 3  maxs ~= +Inf ?
            elseif mins(ix,jx)==0 && maxs(ix,jx)>0 && maxs(ix,jx)<10000
                couplings(ix,jx) = 3;
            % directionally (j->i) coupled - 4 
            elseif mins(ix,jx)>0 && maxs(ix,jx) > 10000
                couplings(ix,jx) = 4;
            end
        end
    end
end
end