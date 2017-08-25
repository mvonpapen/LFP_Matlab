Projects: 
 - Analyses of apomorphine effects (involving phase coherence classification (PCC))
 - Analyses of tremor patients (involving Granger causality (GC))
 - Resting state network (RSN) analyses of M/EEG (data example in ./matlab/brainstorm3/rs_BayHan1_v0.mat)

Folder structure:
  - ARCHIVE: results and scripts that were deemed to be archive-worthy, incl. a backup_file
  - DATA: data
  - cheops: scripts and data taken from cheops, GC of tremor and restng state network analyses
  - matlab_afs: scripts from neuro-project's afs folder. these scripts were tested and work well
  - matlab: large folder containing everything
    - AnalysisApoPatients: pics and results for some patients
    - brainstorm3: brainstorm related content for RSN AnalysisApoPatients
    - Brown2015: scripts to show that Brown2015 results may only show noise correlation
    - Grangwer: scripts and results of GC tremor patient analysis
      - Granger_Detteo: scripts from Detteo et al publication
      - pics: results
      - Tutorial: tutorial from Benjamin Mille
    - moved_to_afs: functions that have been taken out of matlab path, because they have been copied to matlab_afs
    - RSN: scripts and results of RSN analyses
    - shadedErrorBar: shadedErrorBar function
    - temp: temporary stuff
    - Test_pahe_thresh: scripts and data used to calculate phase/time/freq. resolution
    - toolbox: scripts from specific toolboxes that were at the moment not available and thus copied from web

Most important scripts and mat-files:
  - calc_coherence_par.m: calculates coherence from data files (in matlab)
  - calc_pow_v2.m: invoked by former, does the real calculation (in matlab)
  - data_[task]_v3.mat: latest data files (in matlab)
  - load_dat_from_DATA.m: loads data from raw files (in matalab_afs)


Authors: Michael von Papen, Felix Gerick, Benjamin Mille
Date: 02.08.2017