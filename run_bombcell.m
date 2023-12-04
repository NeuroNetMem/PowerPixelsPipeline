function run_bombcell(ephysKilosortPath, ephysBinFile, ephysMetaFile, savePath, decompressDataLocal)

 
    %% set paths - EDIT THESE 
    %ephysKilosortPath = 'D:\NeuropixelData\450409\20231012\raw_ephys_data\probe00\kilosort2_5\sorter_output';% path to your kilosort output files 
    %ephysRawDir = dir('D:\NeuropixelData\450409\20231012\raw_ephys_data\probe00\_spikeGLX_ephysData_g2_t0.imec0.ap.bin'); % path to yourraw .bin or .dat data
    %ephysMetaDir = dir('D:\NeuropixelData\450409\20231012\raw_ephys_data\probe00\_spikeGLX_ephysData_g2_t0.imec0.ap.meta'); % path to your .meta or .oebin meta file
    %saveLocation = 'D:\NeuropixelData\450409\20231012\raw_ephys_data\probe00\kilosort2_5\bombcell_qc'; % where you want to save the quality metrics 
    %savePath = fullfile(saveLocation, 'qMetrics'); 
    %decompressDataLocal = 'D:\NeuropixelData\450409\20231012\raw_ephys_data\probe00\'; % where to save raw decompressed ephys data 
    

    %% load data 
    [spikeTimes_samples, spikeTemplates, templateWaveforms, templateAmplitudes, pcFeatures, ...
        pcFeatureIdx, channelPositions] = bc_loadEphysData(ephysKilosortPath);
    
    %% detect whether data is compressed, decompress locally if necessary
    rawFile = bc_manageDataCompression(dir(ephysBinFile), decompressDataLocal);
    
    %% which quality metric parameters to extract and thresholds 
    param = bc_qualityParamValues(dir(ephysMetaFile), rawFile, ephysKilosortPath); 
    param.unitType_for_phy = 1;
    param.plotGlobal = 0;
    % param = bc_qualityParamValuesForUnitMatch(ephysMetaDir, rawFile) % Run this if you want to use UnitMatch after
    
    %% compute quality metrics 
    rerun = 0;
    qMetricsExist = ~isempty(dir(fullfile(savePath, 'qMetric*.mat'))) || ~isempty(dir(fullfile(savePath, 'templates._bc_qMetrics.parquet')));
    
    if qMetricsExist == 0 || rerun
        [~, ~] = bc_runAllQualityMetrics(param, spikeTimes_samples, spikeTemplates, ...
            templateWaveforms, templateAmplitudes,pcFeatures,pcFeatureIdx,channelPositions, savePath);
    else
        [~, ~] = bc_loadSavedMetrics(savePath); 
    end
end
