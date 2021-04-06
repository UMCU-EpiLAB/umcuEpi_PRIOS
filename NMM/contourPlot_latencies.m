% Contour plot
function contourPlot_latencies(NMM, myDataPath)

mode = {'N1-latency', 'P1-latency','N2-latency'};
 

for J = 1:size(mode,2)
    
    figure('Position',[7,7,948,1055]);
    
   if strcmp(mode{J},'N1-latency')
       x_cor = NMM.outcome.N1_lat;                             % saved after adapt_gain.m
       lat_peak = 28.0;
       xval = abs((x_cor/lat_peak)-1);               % xval with values of the desirec N1-latency will be 1, the further the deviation, the more different colour.
       
   elseif strcmp(mode{J},'P1-latency')
       x_cor = NMM.outcome.P1_lat;
       lat_peak =  65.7;
       xval = abs((x_cor/lat_peak)-1);               
       
   elseif strcmp(mode{J},'N2-latency')
       x_cor = NMM.outcome.N2_lat;
       lat_peak = 219.7;
       xval = abs((x_cor/lat_peak)-1);               

   end
   
    contourf(xval)
    hold on
    % Plot latency value in the plot.
%     N = size(FI_reci,2);
%     x = repmat(1:N,N,1);                                          % generate x-coordinates
%     y = x';                                                       % generate y-coordinates

%     t = num2cell(round(x_cor));                                   % extact values into cells
%     t = cellfun(@num2str, t, 'UniformOutput', false);             % convert to string  
%     text(x(:), y(:), t, 'HorizontalAlignment', 'Center','FontSize',9)

    colormap(autumn)
    ax = gca;
    set(gca,'YDir','normal','ColorScale','log')

    ax.XTick = [1:2:size(x_cor,2)] ;
    ax.XTickLabel = NMM.settings.FI_timeC(1:2:end);
    ax.FontSize = 12;
    ax.YTick =[1:2:size(x_cor,2)] ;
    ax.YTickLabel = NMM.settings.FI_gain(1:2:end);

    xlabel('Time constant of Fast Inhibitory population','FontSize',12,'FontWeight','bold')
    ylabel('Gain of Fast Inhibitory population','FontSize',12,'FontWeight','bold')
    cbh = colorbar;

    title(sprintf('Simulations %s propofol-SPES (in vivo propofol-SPES = %1.1f ms)',mode{J},lat_peak))

    % Set label for colorbar
    new_bar_label = char();
    for i = 1:size(cbh.Ticks,2)
        new_bar_label(i,:) = sprintf('%1.2f \x00B1 %1.3f',lat_peak,(cbh.Ticks(i)*lat_peak));
    end

    set(cbh,'TickLabels', new_bar_label)



     % Save figures
    outlabel=sprintf('ContourPlot_%s.png',mode{J});
    path = fullfile(myDataPath.save_fig_loc);
    
    if ~exist(path, 'dir')
        mkdir(path);
    end
    saveas(gcf,[path,outlabel],'png')
    
end
end

