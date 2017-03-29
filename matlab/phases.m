%% find phase angles from cross-wavelet transformation Wxy in a
%% specified frequency band using only coherent intervals

function phase = phases ( f, Wxy, frange, varargin )

args=struct('coi', 0, ...
            'sig', 0.5, ...  %p05=[0.86 0.645 0.5] for tsmo=[1 2 3]
            'Cxy', 0, ...
            'avg', 0 );
        
args=parseArgs(varargin,args);

if args.coi~=0
    Wxy=coi2nan(f,Wxy,args.coi); 
end

fi = f>=min(frange) & f<=max(frange);
Wxy = Wxy(fi,:);
if args.avg==0
    Wxy = Wxy(:);
end

if args.Cxy~=0
    if args.coi~=0
        Cxy=coi2nan(f,args.Cxy,args.coi);
    end
    Cxy=Cxy(fi,:);
    Cxy=Cxy(:);
    Wxy=Wxy(Cxy>=args.sig);
end

phase=angle(Wxy);

if args.avg==1
    phase=nanmean(phase,1);
end