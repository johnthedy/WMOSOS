clc;clear;tic;close all

%% WMOSOS optimization

% This is multiobjective optimization algorithm modified from SOS
% Current algorithm can only be applied into two objectives problem
% increment = objective weighting interval
% ratiolist = objective weighting assign into each organism
% ite = total iteration number
% ecosize = organism number depends on weighting interval
% lb = lower bound of each optimized parameters
% ub = upper bound of each optimized parameters
% fitness = consist of fobj 1 and fobj 2 (first and second objective) of each eco
% fitnessr = combined fitness of each eco after multiplied by weighting in ratiolist
% PF = Pareto Front Objective
% PFeco = eco corresponding to PF

%% Input
increment=0.05;
ite=20;
n=2;
lb=[-10 -10];
ub=[10 10];

%% START WMOSOS
nd=0:increment:1;
ratiolist=[nd' 1-(nd')];
ecosize=length(ratiolist);
eco=zeros(ecosize,n);
fitnessr=zeros(ecosize,1);
fitness=zeros(ecosize,2);
for i=1:ecosize
    ratio=ratiolist(i,:);
    eco(i,:)=rand(1,n).*(ub-lb)+lb;
    f=fobj(eco(i,:));
    PF(i,:)=f;
    fitnessr(i,:)=sum(f.*ratio);
    fitness(i,:)=f;
end

DOMINATED  = checkDomination(PF);
PF=PF(~DOMINATED,:);
PFeco=eco(~DOMINATED,:);

for h=1:ite
    disp(h)
    for i=1:ecosize

        % Update the best Organism
        r=1/(1+exp(-(10*(h-1)/(0.5*ite-1)-10)));
        if r<0.5
            Usedfitness=fitness;
            Usedfitnesseco=eco;
        else
            Usedfitness=PF;
            Usedfitnesseco=PFeco;
        end

        [szA,~]=size(Usedfitness);
        front=1:1:szA;
        distance = crowdingDistance(front,Usedfitness);
        [d1,d2]=sort(distance,'descend');
        d3=d2(~isinf(d1));
        if isempty(d3)==1
            bestOrganism=Usedfitnesseco(1,:);
        else
            bestOrganism=Usedfitnesseco(d3(1),:);
        end
        
        if length(PF)>ecosize
            [szA,~]=size(PF);
            front=1:1:szA;
            distance = crowdingDistance(front,PF);
            [d1,d2]=sort(distance,'descend');
            PF=PF(d2(1:ecosize),:);
            PFeco=PFeco(d2(1:ecosize),:);
        end
    
        if rand()<0.5
            %Mutualism Phase
            j=i;
            while i==j
                seed=randperm(ecosize);
                j=seed(1);
            end
            % Determine Mutual Vector & Beneficial Factor
            mutualVector=mean([eco(i,:);eco(j,:)]);
            BF1=round(1+rand); BF2=round(1+rand);
            % Calculate new solution after Mutualism Phase
            ecoNew1=eco(i,:)+rand(1,n).*(bestOrganism-BF1.*mutualVector);
            ecoNew1=bound(ecoNew1,ub,lb);
            % Evaluate the fitness of the new solution
            f=fobj(ecoNew1);
            PF=vertcat(PF,f);
            PFeco=vertcat(PFeco,ecoNew1);
            DOMINATED  = checkDomination(PF);
            PF=PF(~DOMINATED,:);
            PFeco=PFeco(~DOMINATED,:);

            for k=1:ecosize
                ratio=ratiolist(k,:);
                y=sum(f.*ratio);
                if y<fitnessr(k)
                    fitnessr(k)=y;
                    fitness(k,:)=f;
                    eco(k,:)=ecoNew1;
                end
            end
            % End of Mutualism Phase
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            % Commensialism Phase
            j=i;
            while i==j
                seed=randperm(ecosize);
                j=seed(1);
            end
            % Calculate new solution after Commensalism Phase
            ecoNew1=eco(i,:)+(rand(1,n)*2-1).*(bestOrganism-eco(j,:));
            ecoNew1=bound(ecoNew1,ub,lb);
            % Evaluate the fitness of the new solution
            f=fobj(ecoNew1);
            PF=vertcat(PF,f);
            PFeco=vertcat(PFeco,ecoNew1);
            DOMINATED  = checkDomination(PF);
            PF=PF(~DOMINATED,:);
            PFeco=PFeco(~DOMINATED,:);

            % Accept the new solution if the fitness is better
            for k=1:ecosize
                ratio=ratiolist(k,:);
                y=sum(f.*ratio);
                if y<fitnessr(k)
                    fitnessr(k)=y;
                    fitness(k,:)=f;
                    eco(k,:)=ecoNew1;
                end
            end
            % End of Commensalism Phase
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Parasitism Phase
        j=i;
        while i==j
            seed=randperm(ecosize);
            j=seed(1);
        end
        % Determine Parasite Vector & Calculate the fitness
        parasiteVector=eco(i,:);
        seed=randperm(n);
        pick=seed(1:ceil(rand*n));  % select random dimension
        parasiteVector(:,pick)=rand(1,length(pick)).*(ub(pick)-lb(pick))+lb(pick);
        f=fobj(parasiteVector);
        PF=vertcat(PF,f);
        PFeco=vertcat(PFeco,parasiteVector);
        DOMINATED  = checkDomination(PF);
        PF=PF(~DOMINATED,:);
        PFeco=PFeco(~DOMINATED,:);

        % Kill organism j and replace it with the parasite
        % if the fitness is lower than the parasite
        for k=1:ecosize
            ratio=ratiolist(k,:);
            y=sum(f.*ratio);
            if y<fitnessr(k)
                fitnessr(k)=y;
                fitness(k,:)=f;
                eco(k,:)=parasiteVector;
            end
        end
        % End of Parasitism Phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    clf
    hold on
    scatter(PF(:,1),PF(:,2),'filled')
    grid on
    hold off
    xlabel('fobj1')
    ylabel('fobj2')
    pause(0.1)
end




