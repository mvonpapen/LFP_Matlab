%%  Save time series for given LFP or EMG channel with additional information
%%  in header
%%
%%  Input:
%%      t       time
%%      data    signal
%%
%%  Optional:
%%      type    LFP or EMG
%%      fname   file name (EMGxxx)
%%      pname   name of patient
%%      site    site number
%%      T       minmax(t) of data
%%      fhead   string describing the filter method
%%      chstr   string describing the channel (C,M,A,L,P,EDC,FDI,FDL)
%%      dt      Time increment
%%      rect    1 if data is rectified, 0 if not

function fid = save_ts ( t, data, varargin )


args=struct('type','n/a',...
            'fname','n/a',...
            'pname','n/a',...
            'site',0,...
            'T',[0 0]
            'fhead','n/a',...
            'chstr','n/a',...
            'dt',1/2500,...
            'rect',0);
        
args=parseArgs(varargin,args);
