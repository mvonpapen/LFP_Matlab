%% Compute and plot averaged spectra
function [md1, md2] = avg_spectra2 ( fnam, datastr1, datastr2, varargin )

args = struct('noplot', false, ...
              'pat',    [], ...
              'fit',    [10 90], ...
              'log',    0, ...
              'data',   []);
args   = parseArgs(varargin, args);

fit    = args.fit; % fit range in Hz
pat    = args.pat;
noplot = args.noplot;
log    = args.log;
data   = args.data;
lab    = {'OFF', 'ON'};


% Load data
if isempty(data)
    data  = load(fnam);
    data1 = data.(datastr1);
    data2 = data.(datastr2);
    f     = data.f;
else
    f     = data{1};
    data1 = data{2};
    data2 = data{3};
end


if ~isempty(pat)
    if ~iscell(pat)
        data1 = data1(:,:,pat);
        data2 = data2(:,:,pat);
    else
        data1 = data1(:,:,pat{1});
        data2 = data2(:,:,pat{2});
    end
end

if log==1
    data1(data1==0) = 1e-10;
    data2(data2==0) = 1e-10;
    data1 = log10(data1);
    data2 = log10(data2);
end


% Parameter
fi             = find(f>=fit(1) & f<=fit(2) & ~(f>45 & f<55));
[~,  nc1, np1] = size(data1);
[nf, nc2, np2] = size(data2);


% Compute average spectra for each data set
md1 = squeeze(nanmean(data1(fi,:,:)));
md2 = squeeze(nanmean(data2(fi,:,:)));
md1 = permute(repmat(md1,1,1,nf),[3 1 2]);
md2 = permute(repmat(md2,1,1,nf),[3 1 2]);
if log~=1
    md1 = data1(:,:)./md1(:,:);
    md2 = data2(:,:)./md2(:,:);
else
    md1 = data1(:,:)-md1(:,:);
    md2 = data2(:,:)-md2(:,:);
end


% Plot results
N1 = sum(~isnan(md1),2)';
N2 = sum(~isnan(md2),2)';
lineprops.col={'r','g'};
if ~noplot
    if log~=1
        mseb( f, [nanmean(md1,2)'; nanmean(md2,2)'], ...
            [nanstd(md1,[],2)'./sqrt(N1); nanstd(md2,[],2)'./sqrt(N2)], lineprops)
%         h1 = shadedErrorBar(f,nanmean(md1,2), ...
%             [nanstd(md1,0,2)'; min([nanstd(md1,0,2)'./sqrt(N1'); nanmean(md1,2)'-1e-20])], 'r', 1);
        set(gca, 'YSca', 'log');
%         hold all
%         h2 = shadedErrorBar(f,nanmean(md2,2), ...
%             [nanstd(md2,0,2)'; min([nanstd(md2,0,2)'./sqrt(N2'); nanmean(md2,2)'-1e-20])], 'g', 1);
    else
%         h1 = shadedErrorBar(f,nanmean(md1,2), nanstd(md1,0,2)./sqrt(N1), 'r', 1);
%         hold all
%         h2 = shadedErrorBar(f,nanmean(md2,2), nanstd(md2,0,2)./sqrt(N2), 'g', 1);
        mseb( f, [nanmean(md1,2)'; nanmean(md2,2)'], ...
            [nanstd(md1,[],2)'./sqrt(N1); nanstd(md2,[],2)'./sqrt(N2)], lineprops)
    end
%     legend([h1.patch h2.patch], lab{1}, lab{2})
    legend(lab{1}, lab{2})
end