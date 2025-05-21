
function [sp_net_community,sp_net_leftover,sp_net_ind,NetInteraction,Quartet,result_one]=quartetIdent(speciesList,common_richness)


% initCobraToolbox;
% changeCobraSolver('ibm_cplex');

% Calculate the growth rate for varying richness
[microbe_alone_GR, richness] = calcMicrobeAloneGr(speciesList, common_richness);

% Find the hidden cooperation
[result_one]=Hid_cooperation(speciesList, common_richness);

% Identifying the left-over abundances
if size(result_one.BM,1)~=0 && size(result_one.BM,2)~=0
    abun_c(:,1)=result_one.vBM;
    [~,index_rich_bhu]=min(abs((richness-common_richness)));
    abun_ind(:,1)=microbe_alone_GR(:,index_rich_bhu);
    
    richness_ind(:,1)=sum_ind_rich(abun_c,microbe_alone_GR,richness);
    
    for inIterNo1 = 1:2
        
        left_over_env(inIterNo1,1)=(common_richness-(sum(richness_ind)-richness_ind(inIterNo1)));
        if left_over_env(inIterNo1,1)<0
            left_over_env(inIterNo1,1)=0;
        end
        
        [~,index_rich(inIterNo1,1)]=min(abs((richness-left_over_env(inIterNo1,1))));
        
        if index_rich(inIterNo1,1)==1
            abun_left_over(inIterNo1,1)=1*0;
        else
            abun_left_over(inIterNo1,1)=microbe_alone_GR(inIterNo1,index_rich(inIterNo1,1));
        end
        if abun_c(inIterNo1,1)~=0
            sp_net_community(inIterNo1)=(abun_c(inIterNo1,1));
            sp_net_ind(inIterNo1)=(abun_ind(inIterNo1,1));
            sp_net_leftover(inIterNo1)=(abun_left_over(inIterNo1,1));
        else
            sp_net_community(inIterNo1)=0;
            sp_net_ind(inIterNo1)=0;
            sp_net_leftover(inIterNo1)=0;
        end
        
    end
    
    netAB=sp_net_community(1)-sp_net_ind(1);
    netAB_pos=sp_net_community(1)-sp_net_leftover(1);
    netAB_neg=netAB-netAB_pos;
    
    netBA=sp_net_community(2)-sp_net_ind(2);
    netBA_pos=sp_net_community(2)-sp_net_leftover(2);
    netBA_neg=netBA-netBA_pos;
    
    NetInteraction = [netAB, netBA];
    Quartet(1)=netAB_pos;
    Quartet(2)=netAB_neg;
    Quartet(3)=netBA_pos;
    Quartet(4)=netBA_neg;
end
end



function [result_orig] = Hid_cooperation(speciesList, common_richness)

numOfSp = length(speciesList);
foldIncrease=ones(1,numOfSp)*100*1;

% Find the biomass indices of the species
Sp_bio_ind = zeros(numOfSp,1);
Sp_bio_ind(1,1) = find(speciesList{1}.c==1);
Sp_bio_ind(2,1) = find(speciesList{2}.c==1);


% Find the ATP maintenance reaction ID of the species
Sp_atp_maint = zeros(numOfSp,1);
for iterNo = 1:numOfSp
    if findRxnIDs(speciesList{iterNo},'DM_atp_c_')
        Sp_atp_maint(iterNo,1) = findRxnIDs(speciesList{iterNo},'DM_atp_c_');
    elseif findRxnIDs(speciesList{iterNo},'ATPM')
        Sp_atp_maint(iterNo,1) =  findRxnIDs(speciesList{iterNo},'ATPM');
    end
end

% Combine the models in a community
modelJoint1 = createMultipleSpeciesModel(speciesList);

% Find the number of cross-feeding reactions and the internal reactions

tot_rxns = length(speciesList{1}.rxns) + length(speciesList{2}.rxns);

tot_exrxns = length(modelJoint1.rxns)-tot_rxns;

% Set the exchange reaction lower bounds
for exInd = 1:tot_exrxns
    modelJoint1.lb(exInd+tot_rxns)=0;
    for iterNo = 1:numOfSp
        if  findRxnIDs(speciesList{iterNo},strrep(modelJoint1.rxns(exInd+tot_rxns),'[u]','(e)'))~=0
            modelJoint1.lb(exInd+tot_rxns)=speciesList{iterNo}.lb(findRxnIDs(speciesList{iterNo},strrep(modelJoint1.rxns(exInd+tot_rxns),'[u]','(e)')));
        elseif findRxnIDs(speciesList{iterNo},strrep(modelJoint1.rxns(exInd+tot_rxns),'[u]','_e'))~=0
            modelJoint1.lb(exInd+tot_rxns)=speciesList{iterNo}.lb(findRxnIDs(speciesList{iterNo},strrep(modelJoint1.rxns(exInd+tot_rxns),'[u]','_e')));

        end
    end
end

%Making the medium rich based on the user input
for exInd = 1:tot_exrxns
    modelJoint1.lb((exInd+tot_rxns))=modelJoint1.lb((exInd+tot_rxns))*common_richness;
    modelJoint1.ub((exInd+tot_rxns))=modelJoint1.ub((exInd+tot_rxns))*common_richness;
end


%transferring individual constraints
for modelNo=1:numOfSp
    model=speciesList{modelNo,1};
    stringg1=strcat('model',num2str(modelNo));
    stringg1=strcat(stringg1,'_IEX');
    
    for rxnNo=1:length(model.rxns)
        reaction=model.rxns(rxnNo);
        if isempty(strfind(reaction{1,1},'EX_'))==0 && isempty(strfind(reaction{1,1},'EX_biomass'))==1 && isempty(strfind(reaction{1,1},'EX_adpcbl(e)'))==1
            stringg2=strrep(model.rxns(rxnNo),'EX',stringg1);
            stringg3=strrep(stringg2{1,1},'(e)','[u]tr');
            if findRxnIDs(modelJoint1,stringg3)>0
                modelJoint1.lb(findRxnIDs(modelJoint1,stringg3))=model.lb(rxnNo)*foldIncrease(modelNo);
            end
        else
            
        end
    end
end

% Set the biomass index
string_gen_biomass = cell(numOfSp,1);
string_gen_maintenance = cell(numOfSp,1);
for iterNo = 1:numOfSp
    model=speciesList{iterNo,1};
    %      string_gen_biomass{iterNo,1}=strcat('model',);
    %     string_gen_biomass{iterNo,1}=strcat('model',num2str(iterNo),'_');
    string_gen_biomass{iterNo,1}=strcat('model',num2str(iterNo));
    string_gen_biomass{iterNo,1}=strcat(string_gen_biomass{iterNo,1},'_');
    string_gen_maintenance{iterNo,1}='ATPM';
    if findRxnIDs(modelJoint1,strcat(string_gen_biomass{iterNo,1},model.rxns(Sp_bio_ind(iterNo))))~=0
        modelJoint1.c(findRxnIDs(modelJoint1,strcat(string_gen_biomass{iterNo,1},model.rxns(Sp_bio_ind(iterNo)))))=1;
    end
end


% Setting the variables for SteadyCom
modelJoint1.csense = char('E' * ones(1,numel(modelJoint1.mets))); % correct the csense
[modelJoint1.infoCom, modelJoint1.indCom] = getMultiSpeciesModelId(modelJoint1, string_gen_biomass);

index=find(modelJoint1.c==1);
options.spBm =modelJoint1.rxns(index);
options.spAbbr = string_gen_biomass;
options.spATPM = string_gen_maintenance;
options.metExId = '[e]';
modelJoint1.infoCom.spBm = modelJoint1.rxns(index); % .spBm for organism biomass reactions
modelJoint1.indCom.spBm = index;

% Run SteadyCom to obtain the community biomass
[~, result_orig, ~, ~, ~] = SteadyCom(modelJoint1);

end

function [microbe_alone_GR, richness] = calcMicrobeAloneGr(speciesList, common_richness)
% Identify the growth rates of individual organisms for varying richness
% (Growth-richness curve)
numOfSp = length(speciesList);
endVal = 1000;
microbe_alone_GR = zeros(numOfSp,endVal);
richness = zeros(endVal,1);
for iterVal = 1: numOfSp
    model = speciesList{iterVal};
    exRxns = findExcRxns(model);
    exRxnsInd = find(exRxns);
    totExRxns = sum(exRxns);
    for iterNo3= 1:1000
        modelNew = model;
        richness(iterNo3)=3*common_richness*iterNo3/endVal;
        for iterNo4 = 1:totExRxns
            modelNew.lb(exRxnsInd(iterNo4))=modelNew.lb(exRxnsInd(iterNo4))*richness(iterNo3);
        end
        result = optimizeCbModel(modelNew, 'max');
        microbe_alone_GR(iterVal,iterNo3)=result.f;
        
    end
end
end

function richness_individual=sum_ind_rich(abun,microbe_alone_GR,richness)
richness_individual=zeros(size(microbe_alone_GR,1),1);
for numOfSp=1:size(microbe_alone_GR,1)
    
    [~,hb]=min(abs(microbe_alone_GR(numOfSp,:)-abun(numOfSp)).^2);
    %   projected_GR(m)=microbe_alone_GR(hb);
    if richness(hb)~=0.006
        richness_individual(numOfSp)=richness(hb);
    end
end
end