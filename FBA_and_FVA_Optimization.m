%% Initialization of COBRA Toolbox
 
initCobraToolbox(false) % Initialization

changeCobraSolver('gurobi') % Changing solver
%% FBA and FVA optimization

model = readCbModel('Model.xml') % Loading our model
original_model = model           % Preserving original model
weights = readmatrix('M_Protein_weights.xlsx') % Loading in our penalty weights

optima = zeros(1,14)
expressionMatrix = zeros(length(model.rxns),14)

for i = 1:14;
    model = original_model
    currentVector = weights(:,i);
 
    % Setting our biomass reaction bounds
    
    model = changeRxnBounds(model,'Biomass_lf',0.01,'b');
    model = changeRxnBounds(model,'Biomass_st',0.01,'b');
    model = changeRxnBounds(model,'Biomass_rt',0.01,'b');
    model = changeRxnBounds(model,'NT_Biomass_lf',0.01,'b');
    model = changeRxnBounds(model,'NT_Biomass_st',0.01,'b');
    model = changeRxnBounds(model,'NT_Biomass_rt',0.01,'b');
  
    % Setting our maintenance requirements
    
    model = addReaction(model,'ATP_maint_lf','reactionFormula','ATP[cl] -> ADP[cl] + Pi[cl]')
    model = addReaction(model,'ATP_maint_st','reactionFormula','ATP[cs] -> ADP[cs] + Pi[cs]')
    model = addReaction(model,'ATP_maint_rt','reactionFormula','ATP[cr] -> ADP[cr] + Pi[cr]')
    model = addReaction(model,'NT_ATP_maint_lf','reactionFormula','ATP[ncl] -> ADP[ncl] + Pi[ncl]')
    model = addReaction(model,'NT_ATP_maint_st','reactionFormula','ATP[ncs] -> ADP[ncs] + Pi[ncs]')
    model = addReaction(model,'NT_ATP_maint_rt','reactionFormula','ATP[ncr] -> ADP[ncr] + Pi[ncr]')
    model = changeRxnBounds(model,'ATP_maint_lf',0.03635,'b')
    model = changeRxnBounds(model,'ATP_maint_st',0.03635,'b')
    model = changeRxnBounds(model,'ATP_maint_rt',0.03635,'b')
    model = changeRxnBounds(model,'NT_ATP_maint_lf',0.03635,'b')
    model = changeRxnBounds(model,'NT_ATP_maint_st',0.03635,'b')
    model = changeRxnBounds(model,'NT_ATP_maint_rt',0.03635,'b')
    model = addReaction(model,'NADPH_maint_lf','reactionFormula','NADPH[cl] -> NADP[cl] + H[cl]')
    model = addReaction(model,'NADPH_maint_st','reactionFormula','NADPH[cs] -> NADP[cs] + H[cs]')
    model = addReaction(model,'NADPH_maint_rt','reactionFormula','NADPH[cr] -> NADP[cr] + H[cr]')
    model = addReaction(model,'NT_NADPH_maint_lf','reactionFormula','NADPH[ncl] -> NADP[ncl] + H[ncl]')
    model = addReaction(model,'NT_NADPH_maint_st','reactionFormula','NADPH[ncs] -> NADP[ncs] + H[ncs]')
    model = addReaction(model,'NT_NADPH_maint_rt','reactionFormula','NADPH[ncr] -> NADP[ncr] + H[ncr]')
    model = changeRxnBounds(model,'NADPH_maint_lf',0.0128','b')
    model = changeRxnBounds(model,'NADPH_maint_st',0.0128','b')
    model = changeRxnBounds(model,'NADPH_maint_rt',0.0128','b')
    model = changeRxnBounds(model,'NT_NADPH_maint_lf',0.0128','b')
    model = changeRxnBounds(model,'NT_NADPH_maint_st',0.0128','b')
    model = changeRxnBounds(model,'NT_NADPH_maint_rt',0.0128','b')

    % Setting our inter-tissue transport ratios
    
    model = addRatioReaction(model,{'Ext_Ala_c_lf_cpls','NT_Ext_Ala_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Arg_c_lf_cpls','NT_Ext_Arg_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Asn_c_lf_cpls','NT_Ext_Asn_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Asp_c_lf_cpls','NT_Ext_Asp_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Cys_c_lf_cpls','NT_Ext_Cys_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Gln_c_lf_cpls','NT_Ext_Gln_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Glu_c_lf_cpls','NT_Ext_Glu_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Gly_c_lf_cpls','NT_Ext_Gly_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_His_c_lf_cpls','NT_Ext_His_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Ile_c_lf_cpls','NT_Ext_Ile_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Leu_c_lf_cpls','NT_Ext_Leu_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Lys_c_lf_cpls','NT_Ext_Lys_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Met_c_lf_cpls','NT_Ext_Met_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Phe_c_lf_cpls','NT_Ext_Phe_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Pro_c_lf_cpls','NT_Ext_Pro_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Ser_c_lf_cpls','NT_Ext_Ser_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Thr_c_lf_cpls','NT_Ext_Thr_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Trp_c_lf_cpls','NT_Ext_Trp_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Tyr_c_lf_cpls','NT_Ext_Tyr_c_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Val_c_lf_cpls','NT_Ext_Val_c_lf_cpls'},[1 3]);

    model = addRatioReaction(model,{'Ext_Ala_c_lf_cpsr','NT_Ext_Ala_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Arg_c_lf_cpsr','NT_Ext_Arg_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Asn_c_lf_cpsr','NT_Ext_Asn_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Asp_c_lf_cpsr','NT_Ext_Asp_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Cys_c_lf_cpsr','NT_Ext_Cys_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Gln_c_lf_cpsr','NT_Ext_Gln_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Glu_c_lf_cpsr','NT_Ext_Glu_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Gly_c_lf_cpsr','NT_Ext_Gly_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_His_c_lf_cpsr','NT_Ext_His_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Ile_c_lf_cpsr','NT_Ext_Ile_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Leu_c_lf_cpsr','NT_Ext_Leu_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Lys_c_lf_cpsr','NT_Ext_Lys_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Met_c_lf_cpsr','NT_Ext_Met_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Phe_c_lf_cpsr','NT_Ext_Phe_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Pro_c_lf_cpsr','NT_Ext_Pro_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Ser_c_lf_cpsr','NT_Ext_Ser_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Thr_c_lf_cpsr','NT_Ext_Thr_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Trp_c_lf_cpsr','NT_Ext_Trp_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Tyr_c_lf_cpsr','NT_Ext_Tyr_c_lf_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Val_c_lf_cpsr','NT_Ext_Val_c_lf_cpsr'},[1 3]);    
 
    model = addRatioReaction(model,{'Ext_Ala_c_st_cpsr','NT_Ext_Ala_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Arg_c_st_cpsr','NT_Ext_Arg_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Asn_c_st_cpsr','NT_Ext_Asn_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Asp_c_st_cpsr','NT_Ext_Asp_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Cys_c_st_cpsr','NT_Ext_Cys_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Gln_c_st_cpsr','NT_Ext_Gln_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Glu_c_st_cpsr','NT_Ext_Glu_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Gly_c_st_cpsr','NT_Ext_Gly_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_His_c_st_cpsr','NT_Ext_His_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Ile_c_st_cpsr','NT_Ext_Ile_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Leu_c_st_cpsr','NT_Ext_Leu_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Lys_c_st_cpsr','NT_Ext_Lys_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Met_c_st_cpsr','NT_Ext_Met_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Phe_c_st_cpsr','NT_Ext_Phe_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Pro_c_st_cpsr','NT_Ext_Pro_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Ser_c_st_cpsr','NT_Ext_Ser_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Thr_c_st_cpsr','NT_Ext_Thr_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Trp_c_st_cpsr','NT_Ext_Trp_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Tyr_c_st_cpsr','NT_Ext_Tyr_c_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Val_c_st_cpsr','NT_Ext_Val_c_st_cpsr'},[1 3]);
   
    model = addRatioReaction(model,{'Ext_Suc_lf_cpls','NT_Ext_Suc_lf_cpls'},[1 3]);
    model = addRatioReaction(model,{'Ext_Suc_st_cpsr','NT_Ext_Suc_st_cpsr'},[1 3]);
    model = addRatioReaction(model,{'Ext_Suc_lf_cpsr','NT_Ext_Suc_lf_cpsr'},[1 3]);
    
    % For low light condition
    
    model = addRatioReaction(model,{'RBC_h_lf', 'RBO_h_lf'},[1 3.5]);                        
    model = addRatioReaction(model,{'DM_DummySucrose_c_lf','DM_DummyStarch_h_lf'},[1 3.8]);   
    
    % For high light condition
    
    %model = addRatioReaction(model,{'RBC_h_lf', 'RBO_h_lf'},[1 2.3]);
    %model = addRatioReaction(model,{'DM_DummySucrose_c_lf','DM_DummyStarch_h_lf'},[1 6.2]);
    
    model = addRatioReaction(model,{'Im_NO3_rt','NT_Im_NO3_rt'},[2 3]);
        
    % Setting up our total flux minimization based on our weights
 
    for p = 1:length(model.rxns);
        model.c(p) = 1;
    end;
    
    for p = 1:length(currentVector);
        model.c(p) = currentVector(p);
    end;
    
    model = convertToIrreversible(model);
 
    for n = 1:length(model.rxns);
        name = char(model.rxns(n));
        last = name(end);
        lastTwo = name([end-1,end])
        if last == 'b';
            model.c(n) = model.c(n)*-1;
        elseif lastTwo == '_r';
            model.c(n) = model.c(n)*-1;
        else
            continue
        end
    end
 
    % Optimization
    
    modelSolFinal = optimizeCbModel(model,'min')  
    
    % Saving our results to the expression matrix
    
    for q = 1:length(modelSolFinal.x);
        expressionMatrix(q,i) = modelSolFinal.x(q);
    end
    
    % Setting our optimal value as a constraint for FVA
    
    reactionList = {}
    for y = 1:length(model.rxns);
        reactionList = [reactionList,model.rxns{y}];
    end

    c_values = []
    for x = 1:length(model.c);
        c_values = [c_values,model.c(x)];
    end
    
    c = model.c
    d = modelSolFinal.f
    ineqS = 'L'
    modelConstrained = addCouplingConstraint(model,reactionList,c_values,d,ineqS)
    
    % FVA optimization
    
    reactionsToCheck = {'RBC_h_lf','SBPA_h_lf','SBPase_h_lf','STK_h_lf_f','STK_h_lf_b','FTK_h_lf_f','FTK_h_lf_b','GAPDH1_h_lf','STK_h_lf_f','STK_h_lf_b','FTK_h_lf_f','FTK_h_lf_b','Ru5PE_h_lf_f','Ru5PE_h_lf_b','R5PI_h_lf_f','R5PI_h_lf_b','Ru5PK_h_lf','RBO_h_lf','PGP_h_lf','GlyDH_m_lf','SGAT_p_lf','GCEAK_h_lf','FBPA_h_lf_f','FBPase_h_lf','PGI_h_lf_f','PGI_h_lf_b','PGM_h_lf_f','PGM_h_lf_b','AGPase_h_lf','FBPA_c_lf_f','FBPA_c_lf_b','FBPase_c_lf','PGI_c_lf_f','PGI_c_lf_b','PGM_c_lf_f','PGM_c_lf_b','UGPase_c_lf_f','UGPase_c_lf_b','S6PS_c_lf','Enol_c_lf_b','Enol_c_lf_f','PyrK_c_lf','PyrDH1_m_lf','CitS_m_lf','cACNDHA_m_lf','iCitDHNAD_m_lf','iCitDHNAD_c_lf','iCitDHNADP_c_lf','iCitDHNADP_h_lf','iCitDHNADP_m_lf','MalDH1_m_lf_f','MalDH1_m_lf_b','PEPC2_c_lf_b','PEPC2_c_lf_f','AlaTA_m_lf_f','AlaTA_m_lf_b','AspAT_c_lf_f','AspAT_c_lf_b','AspAT_m_lf_f','AspAT_m_lf_b','AspAT_h_lf_f','AspAT_h_lf_b','AspAT_p_lf_f','AspAT_p_lf_b','GluSNADP_h_lf_f','GluSNADP_h_lf_b','GluDH3NADP_c_lf_b','GluDH3NADP_c_lf_f','GluSNAD_h_lf','GluSFd_h_lf','GluDH1NAD_m_lf','GluDH1NADP_m_lf_f','GluDH1NADP_m_lf_b','GluSeADH_c_lf','GluSeADH_h_lf','GluSeADH_m_lf','GlnS_c_lf','GlnS_h_lf','GlnS_m_lf','ThrS_h_lf','AsnS_c_lf','Tr_TPT1_lf_f','Tr_TPT1_lf_b','Tr_TPT2_lf_f','Tr_TPT2_lf_b','Tr_TPT3_lf_f','Tr_TPT3_lf_b','Im_CO2_lf_b','Im_CO2_lf_f'}
    
    minFluxList = []
    maxFluxList = []
    
    for q=1:length(reactionsToCheck)
        reactionName = reactionsToCheck(q)
        modelConstrained = changeObjective(modelConstrained,reactionName)
        minFlux = optimizeCbModel(modelConstrained,'min')
        maxFlux = optimizeCbModel(modelConstrained,'max')
        minFluxList = [minFluxList,minFlux.f]
        maxFluxList = [maxFluxList,maxFlux.f]
    end
    
    % Saving our FVA results
    
    minFilename = sprintf('Mergner_Protein_LowLight_MinimumFluxes_%d.xls',i)
    maxFilename = sprintf('Mergner_Protein_LowLight_MaximumFluxes_%d.xls',i)
    
    xlswrite(minFilename,minFluxList)
    xlswrite(maxFilename,maxFluxList)
    
end

% Saving our flux results

xlswrite('Mergner_Protein_LowLight_ExpressionMatrix.xls',expressionMatrix)