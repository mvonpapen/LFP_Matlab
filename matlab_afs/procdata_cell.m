function [ data, W, coi, Pw, ff, Pf ] = procdata_cell( dat, art, varargin )
%PROCDATA_CELL Summary of this function goes here
%   Detailed explanation goes here
if nargin<2
    error('need artifact cell array!');
end
args=struct('filter',[],... % for stop band e.g. [45 55; 90 110]
            'order',5,...
            'freq',logspace(0,2.7,30),...
            'rect',0,...
            'dt',1/2500,...
            'inc',0,...
	    'w0', 6);
args=parseArgs(varargin,args);


if ~iscell(dat)
    [ data, W, coi, Pw, ff, Pf ] = ...
            procdata ( dat, 'filter', args.filter, 'order',args.order, ...
            'freq',args.freq,'rect',args.rect,'inc',args.inc,'art',art,'w0',args.w0);
        
else

    data=cell(size(dat));
    W=cell(size(dat));
    coi=cell(size(dat));
    Pw=cell(size(dat));
    ff=cell(size(dat));
    Pf=cell(size(dat));
    if iscell(dat{1})
        for i=1:length(dat)
            data{i}=cell(size(dat{i}));
            W{i}=cell(size(dat{i}));
            coi{i}=cell(size(dat{i}));
            Pw{i}=cell(size(dat{i}));
            ff{i}=cell(size(dat{i}));
            Pf{i}=cell(size(dat{i}));
            for j=1:length(dat{i})
                [ data{i}{j}, W{i}{j}, coi{i}{j}, Pw{i}{j}, ff{i}{j}, Pf{i}{j} ] = ...
                procdata ( dat{i}{j}, 'filter', args.filter, 'order',args.order, ...
                'freq',args.freq,'rect',args.rect,'inc',args.inc,'art',art{i}{j},'w0',args.w0);
            end
        end
    else
        for i=1:length(dat)
            [ data{i}, W{i}, coi{i}, Pw{i}, ff{i}, Pf{i} ] = ...
            procdata ( dat{i}, 'filter', args.filter, 'order',args.order, ...
            'freq',args.freq,'rect',args.rect,'inc',args.inc,'art',art{i},'w0',args.w0);

        end

    end

end
end
