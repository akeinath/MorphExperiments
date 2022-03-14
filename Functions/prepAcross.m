function across = prepAcross()

    across.isPlaceCell = [];

    across.mfr = [];
    across.maps = [];
    across.envs = [];
    across.mds.angles = [];
    across.mds.driftVariability = [];
    across.mds.driftAmount = [];
    across.mds.stress = [];
    
    across.seqAnalysis.maxDecorr = [];
    across.seqAnalysis.transitionPoint = [];
    
    across.drift.numSPart = [];
    across.drift.mfrRankPart = [];
    across.drift.conRankPart = [];
    
    across.RDMs.partitioned = [];
    across.RDMs.whole = [];
    across.RDMs.envs = [];
    
    across.contextPrediction.rawPV.actual = [];
    across.contextPrediction.zeroedPV.actual = [];
    across.contextPrediction.justTracked.actual = [];
    
    across.environmentPrediction.actual = [];
    across.environmentPrediction.actual_similarity = [];
    
    across.contextPrediction.contextErrorDivisions = [];
    across.contextPrediction.contextErrorDivisions_similarity = [];
    
    across.partSim = [];
    across.sim = [];
    across.simXcon = [];
    across.simXseq = [];
    across.class.d2.r2penalty = [];
    across.class.d2.tr2 = [];
    across.class.d2.shuffle_r2penalty = [];
    across.class.d2.shuffle_tr2 = [];
    across.class.d3.r2penalty = [];
    across.class.d3.tr2 = [];
    across.class.d3.pop.r2penalty = [];
    across.class.d3.pop.tr2 = [];
    across.class.d3.shuffle_r2penalty = [];
    across.class.d3.shuffle_tr2 = [];
    
    across.cellAnimalID = [];

    across.coverage = [];
    
    across.newSHC = [];
end