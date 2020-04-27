function plot_ccep_av(dataBase,cfg)
%
% function ccep_plot_av(average_ccep,tt,n1_peak_sample, n1_peak_amplitude,average_ccep_names,...
%     channel_names,good_channels,myDataPath,bids_sub,bids_ses,bids_task,bids_runs,params)
%
% Function plots average CCEPs across conditions per electrode.
%
% input 
%   average_ccep: electrodes X condition (stim pair) X time
%   tt: time
%   n1_peak_sample: electrodes X condition, [] to not show
%   n1_peak_amplitude: [] or electrodes X condition, [] to not show
%   average_ccep_names: ccep condition (stim pair) names
%   channel_names: names of the channels (size electrodes)
%   good_channels: electrode indices to plot
%   myDataPath:
%   bids_sub:
%   bids_ses:
%   bids_ses:
%   bids_task:
%   bids_runs:
%   save_fig: 1 to save, 0 not to save
%
% Dora Hermes, 2020, Multimodal Neuroimaging Lab, Mayo Clinic
% Dorien van Blooijs, 2020, UMC Utrecht

myDataPath.CCEPpath = '/Fridge/users/sifra/derivatives/CCEP/' ;
myDataPath.dataPath = '/Fridge/chronic_ECoG/';

if isempty(dataBase.ccep.n1_peak_sample)
    n1_peak_sample = NaN(size(dataBase.cc_epoch_sorted_avg,1),size(dataBase.cc_epoch_sorted_avg,2));
    n1_peak_amplitude =  NaN(size(dataBase.cc_epoch_sorted_avg,1),size(dataBase.cc_epoch_sorted_avg,2));
else
    n1_peak_sample = dataBase.ccep.n1_peak_sample;
    n1_peak_amplitude = dataBase.ccep.n1_peak_amplitude;
end
    
tt = dataBase.tt;

elnrs_plot = 1:size(dataBase.ch,1);

for ll = 1:length(elnrs_plot)                   % For the number of electrodes
    el_plot = elnrs_plot(ll);
    %figure('Position',[0 0 700 700]),hold on
    figure, hold on
    for kk = 1:length(dataBase.cc_stimchans)    % For the number of stimulation pairs
        this_ccep_plot = squeeze(dataBase.cc_epoch_sorted_avg(el_plot,kk,:));
        this_ccep_plot(tt>-0.010 & tt<0.010) = NaN;
        
        plot(tt,kk*500+zeros(size(tt)),'Color',[.8 .8 .8])
        plot(tt,kk*500+this_ccep_plot)
        if ~isnan(n1_peak_sample(el_plot,kk))
            plot(tt(n1_peak_sample(el_plot,kk)),dataBase.cc_epoch_sorted_avg(el_plot,kk,n1_peak_sample(el_plot,kk))+kk*500,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',2)
        end
    end
    xlim([-.2 1.5])
    ylim([-500,(kk+2)*500])
    set(gca,'YTick',500*(1:length(dataBase.stimpnames)),'YTickLabel',dataBase.stimpnames)
    title([dataBase.ch{el_plot}])
    
    ylabel('stimulated electrodes')
    xlabel('time(s)')
    
    % add amplitude bar
    plot([0.9 0.9],[1000 1500],'k','LineWidth',2)
    text(0.91,1250,['500 ' native2unicode(181,'latin1') 'V'])
    
    if cfg.save_fig==1
        % create folder to save figures
        if ~exist(fullfile(myDataPath.CCEPpath,'av_ccep_figures',dataBase.sub_label,dataBase.ses_label,dataBase.run_label),'dir')
            mkdir(fullfile(myDataPath.CCEPpath,'av_ccep_figures',dataBase.sub_label,dataBase.ses_label,dataBase.run_label));
        end

        % filename
        figureName = fullfile(myDataPath.CCEPpath,'av_ccep_figures',dataBase.sub_label,dataBase.ses_label,dataBase.run_label,...
            [dataBase.sub_label '_' dataBase.ses_label '_' dataBase.run_label '_incomingCCEP_el' dataBase.ch{el_plot}]);
        set(gcf,'PaperPositionMode','auto')
        print('-dpng','-r300',figureName)
        print('-depsc','-r300',figureName)
    else
        pause
    end
    close all
end

end