  util.info("oposSOM is ready to fly! Starting analysis.")
  util.info("Name:", env$preferences$dataset.name)

  #### Preparation & Calculation part ####
  env <- pipeline.checkInputParameters(env)
  if (!env$passedInputChecking) {
    return()
  }
  
  if(env$preferences$activated.modules$primary.analysis)
  {
    env$preferences$system.info <- Sys.info()
    env$preferences$session.info <- sessionInfo()
    env$preferences$started <- format(Sys.time(), "%a %d %b %Y %X")
  }
  
  if(env$preferences$activated.modules$reporting)
  {
    # create output dirs
    dir.create(paste(env$files.name, "- Results"), showWarnings=FALSE)
    dir.create(paste(env$files.name, "- Results/CSV Sheets"), showWarnings=FALSE)
	setwd(paste(env$files.name, "- Results"))
	
    if(env$preferences$activated.modules$primary.analysis)
    {
      pipeline.qualityCheck(env)
    } 
  }
  if(env$preferences$activated.modules$primary.analysis || env$preferences$activated.modules$geneset.analysis)
  {
    util.info("Loading gene annotation data.")
    env <- pipeline.prepareAnnotation(env)
  }
  
  if(env$preferences$activated.modules$primary.analysis)
  {
    util.info("Processing SOM. This may take several time until next notification.")
    env <- pipeline.prepareIndata(env)
    env <- pipeline.generateSOM(env)
    
    filename <- paste(env$files.name, "pre.RData")
    util.info("Saving environment image:", filename)
    save(env, file=filename)
    
    util.info("Processing Differential Expression Statistics")
    env <- pipeline.diffExpressionStatistics(env)

    util.info("Detecting Spots")
    env <- pipeline.detectSpotsSamples(env)
    env <- pipeline.detectSpotsModules(env)
    env <- pipeline.patAssignment(env)
    env <- pipeline.groupAssignment(env)
  }

  if (env$preferences$activated.modules$geneset.analysis)
  {
    util.info("Calculating Geneset Enrichment")
    env <- pipeline.genesetStatisticSamples(env)
    env <- pipeline.genesetStatisticModules(env)
  }
  
  if (env$preferences$activated.modules$psf.analysis)
  {
    util.info("Calculating Pathway Signal Flow (PSF)")
    env <- pipeline.PSFcalculation(env)    
  }
  
  if(env$preferences$activated.modules$primary.analysis || env$preferences$activated.modules$geneset.analysis)
  {    
    filename <- paste(env$files.name, ".RData", sep="")
    util.info("Saving environment image:", filename)
    save(env, file=filename)
    
    if (file.exists(paste(env$files.name, "pre.RData")) && file.exists(filename))
    {
      file.remove(paste(env$files.name, "pre.RData"))
    }
  }  
    
  #### Reporting part ####
  
  if(env$preferences$activated.modules$reporting)
  {
  
    util.info("Plotting Supporting Information")
    pipeline.supportingMaps(env)
    pipeline.entropyProfiles(env)
    pipeline.topologyProfiles(env)

    
    if(length(env$chromosome.list) > 0)
    {
      util.info("Plotting Chromosome Expression Reports")
      pipeline.chromosomeExpressionReports(env)
    }
    
    if(ncol(env$indata) < 1000)
    {
      util.info("Plotting Sample Portraits")
      pipeline.sampleExpressionPortraits(env)
    } 
    
    if ( env$preferences$activated.modules$sample.similarity.analysis && ncol(env$indata) > 2)
    {    
      util.info("Plotting Sample Similarity Analysis")
      dir.create("Sample Similarity Analysis", showWarnings=FALSE)
      
      pipeline.sampleSimilarityAnalysisED(env)
      pipeline.sampleSimilarityAnalysisCor(env)
      pipeline.sampleSimilarityAnalysisICA(env)
      pipeline.sampleSimilarityAnalysisSOM(env)
    }
    
    if (env$preferences$activated.modules$geneset.analysis)
    {
      dir.create("Geneset Analysis", showWarnings=FALSE)
      
      util.info("Plotting Geneset Enrichment Heatmaps")
      pipeline.genesetOverviews(env)
      
      util.info("Plotting Geneset Profiles and Maps")
      pipeline.genesetProfilesAndMaps(env)
      
      util.info("Calculating Cancer Hallmark Enrichment")
      pipeline.cancerHallmarks(env)
    }
    
    if (env$preferences$activated.modules$psf.analysis)
    {
      util.info("Plotting PSF results")
      pipeline.PSFoutput(env)
    }
    
    util.info("Writing Gene Lists")
    pipeline.geneLists(env)

    util.info("Plotting Summary Sheets (Samples)")
    pipeline.summarySheetsSamples(env)
    
    util.info("Plotting Summary Sheets (Modules & PATs)")
    pipeline.summarySheetsModules(env)
    pipeline.summarySheetsPATs(env)
      

    if(env$preferences$activated.modules$group.analysis && length(unique(env$group.labels)) >= 2)
    {
      util.info("Processing Group-centered Analyses")
      pipeline.groupAnalysis(env)
    }
  
    if(env$preferences$activated.modules$difference.analysis)
    {
      util.info("Processing Difference Analyses")
      pipeline.differenceAnalyses(env)
    }

    util.info("Generating HTML Report")
    pipeline.htmlSampleSummary(env)
    pipeline.htmlModuleSummary(env)
    pipeline.htmlGenesetAnalysis(env)  
    pipeline.htmlPsfAnalysis(env)
    pipeline.htmlSummary(env)
    
  }    
    
  util.info("Finished:", format(Sys.time(), "%a %b %d %X"))
	
	return(env)