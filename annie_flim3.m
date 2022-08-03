% and number of cell export and mode...

clear all

% EXP1
%fld_msk = 'C:\Users\annie\OneDrive - University of Cambridge\Data\Imaging data\FLIM snapshot time course\Clone B1\Rep 1\Segmentation maps\Segmentation maps\';
%fld_fl1 = 'C:\Users\annie\OneDrive - University of Cambridge\Data\Imaging data\FLIM snapshot time course\Clone B1\Rep 1\T1 blanket threshold\T1 blanket threshold\';
% fld_raw = 'E:\Current Users\Annie\4. KRAS mutant cell lines\HPNE\mitoAteam sensor imaging\210603 FLIM snapshot I\210603 FLIM snapshot 1 (raw data)\Data from microscope\';

% EXP2
fld_msk = 'C:\Users\annie\OneDrive - University of Cambridge\Data\Imaging data\FLIM snapshot time course\Clone B1\Rep 2\Segmentation maps\';
fld_fl1 = 'C:\Users\annie\OneDrive - University of Cambridge\Data\Imaging data\FLIM snapshot time course\Clone B1\Rep 2\T1s';

% EXP3
%fld_msk = 'E:\Current Users\Annie\4. KRAS mutant cell lines\HPNE\mitoAteam sensor imaging\210617 FLIM snapshot 3\Ilastik files\Segmentation maps\';
%fld_fl1 = 'E:\Current Users\Annie\4. KRAS mutant cell lines\HPNE\mitoAteam sensor imaging\210617 FLIM snapshot 3\Ilastik files\T1s\';


bNewAnalysis = 0; % %  0 / 1 0: use SPCImage data

mrk_day = {'0d','1d','3d','5d','7d'}; % first sample marker (day)
mrk_smp = {'mc','wt','d','r','v'}; % second sample marker (sample id)

dn = numel(mrk_day);
sn = numel(mrk_smp);

% retrieve file names
file_names = dir([fld_msk '*.tiff']);
file_names = {file_names(~[file_names.isdir]).name}; % purge folder names

fn = numel(file_names);

f1x = 1; f1y = 3;
hf1 = figure;

xhis = (.5:.01:2.1);
hfl1 = zeros(numel(xhis)-1,dn,sn);


sca = {};
scs = {};
for fi=1:fn
    
   % parse file name
   
   base_name     = file_names{fi}(1:strfind(file_names{fi},'_intensity')-1);
   markers       = regexp(base_name,'_','split');
   [dmb day_idx] = ismember(markers{1},mrk_day);
   [dmb smp_idx] = ismember(markers{2},mrk_smp);
   
   % load segmentation
   msk = imread([fld_msk file_names{fi}]); % load mask
   lbl = bwlabel(msk==2);  % label mask
   
   switch bNewAnalysis
       case 0 % load pre analysed data
            fl1 = textread([fld_fl1 base_name '_t1.asc'],'%f'); % load flim1
            fl1 = rot90(reshape(fl1,[256 256]))/1000; % rebuild image and align to segmentation and convert to ns
       case 1 % de novo analysis of the sdt files
            [a b] = read_sdt([fld_raw base_name '.sdt']);
           
   end
   
   on  = max(lbl(:));
    try
        o0 = numel(sca{day_idx}{smp_idx});
    catch
        o0 = 0;
    end
   for oi=1:on
       
       if nnz((lbl==oi).*fl1)>0
       
           sca{day_idx}{smp_idx}(o0+oi) = mean(nonzeros((lbl==oi).*fl1));    % average
           scm{day_idx}{smp_idx}(o0+oi) = median(nonzeros((lbl==oi).*fl1));  % median
           scd{day_idx}{smp_idx}(o0+oi) = mode(nonzeros((lbl==oi).*fl1));  % mode
           scs{day_idx}{smp_idx}(o0+oi) = std(nonzeros((lbl==oi).*fl1));

           his_tmp = histcounts(nonzeros((lbl==oi).*fl1),xhis);
           his_tmp = his_tmp/sum(his_tmp); %each cell will contibute to the final histrogram in equal proportions
           hfl1(:,day_idx,smp_idx,:) = hfl1(:,day_idx,smp_idx,:)+his_tmp' ;
           hfl1(:,day_idx,smp_idx) = hfl1(:,day_idx,smp_idx) / sum(nonzeros(hfl1(:,day_idx,smp_idx)));
       else
           o0=o0-1;
       end
   end
   
  
%    subplot(f1x,f1y,1)
%    imagesc(lbl)
%    axis square 
%    axis off
%    title('segmentation')
% 
%    subplot(f1x,f1y,2)
%    imagesc(fl1)
%    axis square 
%    axis off
%    title('FLIM1')
%    
%    subplot(f1x,f1y,3)
%    boxplot(sca,'PlotStyle','compact')
%    axis square
%    
end

%%
string_analysis = {'means','medians','modes','std'};

for ai=1:3 % do not plot yet std
    switch ai
        case 1
            vals = @(x,y)sca{x}{y};
        case 2
            vals = @(x,y)scm{x}{y};
        case 3
            vals = @(x,y)scd{x}{y};
        case 4
            vals = @(x,y)scs{x}{y};
    end
    %% box plot
    box_values = [];
    box_group  = [];
    for di=1:dn
        for si=1:sn

            box_values = [box_values vals(di,si)];
            box_group  = [box_group (di+(si-1)*dn)*ones(1,numel(vals(di,si)))]

        end
    end


    smp_col = [ repmat([.6 .2 0],[dn 1]); % mCherry
                repmat([ 0  0 0],[dn 1]); % WT
                repmat([ 1  0 0],[dn 1]); % D
                repmat([ 0  1 0],[dn 1]); % R
                repmat([ 0  0 1],[dn 1]); % V
              ];

    figure('name',['box plot - single cell distributions - ' string_analysis{ai}] )
    boxplot(box_values,box_group,'PlotStyle','compact','colors',smp_col)
    set(gca,'xtick',1:25,'xticklabel',{'0','1','3','5','7'})
    xlabel('day post DOX induction')
    ylabel('fluorescence lifetime (ns)')

    text(3,1.7,'mCherry')
    text(8,1.7,'WT')
    text(13,1.7,'G12D')
    text(18,1.7,'G12R')
    text(22,1.7,'G12V')


%% sc distributions

    lt_min = min(box_values);
    lt_max = max(box_values);
    lt_his = (lt_min:0.01:lt_max);
    krn = 6;

    figure('name',['Distribution of single cell ' string_analysis{ai}])
    
    for si=1:sn
        subplot(2,5,si)
        for di=1:dn
            yoffset = (di-1)*.05;
            y = histcounts(vals(di,si),lt_his);
            y = circshift(filter(ones(1,krn)/krn,1,y),-round(krn/2));
            plot((lt_his(1:end-1)+lt_his(2:end))/2,yoffset+y/sum(y),'color',[1-di/dn 0 di/dn .7],'linewidth',3)      
            text(lt_max,yoffset,num2str(numel(vals(di,si))))
            hold on
        end
        title(mrk_smp{si})
        set(gca,'xlim',[lt_min lt_max])
        axis square
    end


    for di=1:dn
        subplot(2,5,di+5)
        for si=1:sn
            yoffset = (si-1)*.08;
            y = histcounts(vals(di,si),lt_his);
            y = circshift(filter(ones(1,krn)/krn,1,y),-round(krn/2));
            plot((lt_his(1:end-1)+lt_his(2:end))/2,yoffset+y/sum(y),'color',[1-si/sn 0 si/sn .7],'linewidth',3)
            text(lt_max,yoffset,num2str(numel(vals(di,si))))
            hold on
        end
        title(mrk_day{di})
        axis square
        set(gca,'xlim',[lt_min lt_max])
    end

end
%% pixel distributions
xval = (xhis(1:end-1)+xhis(2:end))/2;
x1 = 1;
x2 = numel(xval);
krn = 13;


figure('name','pixel distributions')
for si=1:sn
    subplot(2,5,si)
    for di=1:dn
        y = circshift(filter(ones(1,krn)/krn,1,squeeze(hfl1(x1:x2,di,si))),-round(krn/2));
        plot(xval(x1:x2),(di-1)*.035+y,'color',[1-di/dn 0 di/dn .7],'linewidth',3);
        hold on
    end
    axis square
    set(gca,'xlim',[0.9 1.9],'ylim',[-0.01 0.2])
    title(mrk_smp{si})
end

for di=1:dn
    subplot(2,5,di+5)
    for si=1:sn
        y = circshift(filter(ones(1,krn)/krn,1,squeeze(hfl1(x1:x2,di,si))),-round(krn/2));
        ho = plot(xval(x1:x2),(si-1)*.035+y,'color',[1-si/sn 0 si/sn .7],'linewidth',3)
        hold on
    end
    axis square
    set(gca,'xlim',[0.9 1.9],'ylim',[-0.01 0.2])
    title(mrk_day{di})
end
beep


%% sc distributions (standard deviations)

figure('name','Distribution of single cell standard deviations')

lt_min = .003;
lt_max = .33;
lt_his = (lt_min:0.0025:lt_max);
krn = 8;
for si=1:sn
    subplot(2,5,si)
    for di=1:dn
        y = histcounts(scs{di}{si},lt_his);
        y = circshift(filter(ones(1,krn)/krn,1,y),-round(krn/2));
        plot((lt_his(1:end-1)+lt_his(2:end))/2,(di-1)*.03+y/sum(y),'color',[1-di/dn 0 di/dn .7],'linewidth',3)
        hold on
    end
    title(mrk_smp{si})
    set(gca,'xlim',[0 .33])
    axis square
end


for di=1:dn
    subplot(2,5,di+5)
    for si=1:sn
        y = histcounts(scs{di}{si},lt_his);
        y = circshift(filter(ones(1,krn)/krn,1,y),-round(krn/2));
        plot((lt_his(1:end-1)+lt_his(2:end))/2,(si-1)*.04+y/sum(y),'color',[1-si/sn 0 si/sn .7],'linewidth',3)
        hold on
    end
    title(mrk_day{di})
    axis square
    set(gca,'xlim',[0 0.33])
end