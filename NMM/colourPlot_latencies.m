function colourPlot_latencies(NMM, myDataPath)

mode = {'N1-latency', 'P1-latency','N2-latency'};
yAs_gain = NMM.settings.FI_gain;


% Plot the peak-latencies from the NMM 
% The gain is on the y-axis, the variations of the time-constant are
% displayed with various colors. 
% x-axis is the latency (ms)

for J = 1:size(mode,2)
    
    
    figure('Position',[7,7,948,1055]);
    
   if strcmp(mode{J},'N1-latency')
       x_cor = NMM.outcome.N1_lat;                             % saved after adapt_gain.m
       lat_peak = 28.0;
       hold on
       xline(28.0,'-.','LineWidth',3);
       xline(22.7,'-.','LineWidth',3);


   elseif strcmp(mode{J},'P1-latency')
       x_cor = NMM.outcome.P1_lat;
       lat_peak = 65.7;
       hold on
       xline(65.7,'-.','LineWidth',3);
       xline(59.2,'-.','LineWidth',3);
       
   elseif strcmp(mode{J},'N2-latency')
       x_cor = NMM.outcome.N2_lat;
       lat_peak = 219.7;
       hold on
       xline(207.5,'-.','LineWidth',3);
       xline(219.7,'-.','LineWidth',3);

   end

   % plot the result of very latency & gain combination in various colours.
   n=size(x_cor,2);
   c= jet(n);

   for i = 1:size(x_cor,2)
       gain_plot(i,:) = plot(x_cor(:,i),yAs_gain,'Color',c(i,:),'LineWidth',3);
       hold on    
   end

    ylabel('Gain Fast Inhibitory population (mV)','FontSize',12,'FontWeight','bold')
    xlabel(sprintf('%s (ms)',mode{J}),'FontWeight','bold');
    title(sprintf('Simulations %s (in vivo propofol-SPES = %1.1f ms)',mode{J},lat_peak))
        
    ax = gca;
    ax.FontSize = 12;
    ax.YTick = [min(NMM.settings.FI_gain):max(NMM.settings.FI_gain)] ;
    ax.YTickLabel = NMM.settings.FI_gain(1:2:end);
    
      
    % Display legend with value time constant
    g = NMM.settings.FI_timeC;
    legendCell = cellstr(num2str(g', 'g=%-d s^{-1}'));
    grid on
    legend(gain_plot,legendCell,'Location','northeastoutside','Fontsize',12)

    
     % Save figures
    outlabel=sprintf('ColourPlot_1%s.png',mode{J});
    path = fullfile(myDataPath.save_fig_loc);
       
    if ~exist(path, 'dir')
        mkdir(path);
    end
    saveas(gcf,[path,outlabel],'png')
  
    
end

end