function generateOutputFigures(paselected, paquantl, final)

%GENERATEOUTPUTFIGURES Generate output figures - original
%   Prasanna Sritharan, Mario Andrés Muñoz, August 2019
%
% Based on original scripts written by Mario Andrés Muñoz, 2013


% User settings
user = getUserScriptSettings();
outpath = user.OUTPATH1;
limbs = user.LIMBS;

fprintf('Generating output figures.\n');
fprintf('------------------------------------------------\n');

% event timings (from Python OpenSim and SPM workflow)
events = [0, 28.5, 38.25, 68.25, 77, 101];
eventlabels = {'PFO1', 'PFS2', 'NFO1', 'NFS2', 'PFO3', 'PFS4'};
eventlabelspos = [0.5, 28, 38.75, 67.75, 77.50, 100.5];
eventlabelsalign = {'left', 'right', 'left', 'right', 'left', 'right'};

% Limbs
for b=1:2

    fprintf('\nFigures: %s LIMB\n', upper(limbs{b}));    

    % Plot quantiles and PC correlations
    for f=1:length(final.(limbs{b}).labels)
    
        % Get PC info from label
        pcstr = '([^_]+)_(\w+_*\w*)_PC(\d+)';
        toks = regexp(final.(limbs{b}).labels{f},pcstr, 'tokens');
        
        % Get PC index in quantile struct, and info for plot labels
        switch toks{1}{1}
            case 'angle'
                analysis = 'ik';
                ytext1 = 'Angle (deg)';
            case 'moment'
                analysis = 'id';
                ytext1 = 'Moment (%BW*HT)';
        end
        
        fprintf('Generating figure: %s...\n', final.(limbs{b}).labels{f});
        
        figure(f);
        set(gcf, 'Position', [100 100 400 900], 'PaperSize', [29.7 21.0]);
        sgtitle(strrep(final.(limbs{b}).labels{f}, '_', '-'));
        

        % ********************
        % Quantiles
        
        % Data
        data.top.mean = paquantl.(limbs{b}).(analysis).(toks{1}{2}).top.mean(:, str2double(toks{1}{3}));
        data.top.std = paquantl.(limbs{b}).(analysis).(toks{1}{2}).top.std(:, str2double(toks{1}{3}));
        data.top.upper = data.top.mean + data.top.std;
        data.top.lower = data.top.mean - data.top.std;
        data.bottom.mean = paquantl.(limbs{b}).(analysis).(toks{1}{2}).bottom.mean(:, str2double(toks{1}{3}));
        data.bottom.std = paquantl.(limbs{b}).(analysis).(toks{1}{2}).bottom.std(:, str2double(toks{1}{3}));
        data.bottom.upper = data.bottom.mean + data.bottom.std;
        data.bottom.lower = data.bottom.mean - data.bottom.std;
        
        % Plot
        subplot(2, 1, 1);
        hold on;
        qtl = {'top', 'bottom'};
        clrs = {'b', 'r'};
        lsty = {'', '--'};
        for d=1:2
            [ha, hbot, htop] = shadedplot(0:100, data.(qtl{d}).lower', data.(qtl{d}).upper',clrs{d});
            alpha(ha(2), 0.3);
            set(ha, 'HandleVisibility', 'off');   % prevent legend
            delete(hbot);   % remove lower and upper bound lines automatically created by shaded plot
            delete(htop);
            hold on;    % shadedplot sets hold off
            plot(0:100, data.(qtl{d}).mean, [clrs{d} lsty{d}]);        
        end
        for ev=2:5, xline(events(ev), ':'); end
        hold off;
        %legend('25%Q', '75%Q');
        box on;
        xlim([0 100]);
        xlabel('% of step-down-and-pivot task');
        ylim('auto');   % temporary
        ylabel(ytext1);
        yl = ylim;
        for ev=1:6, text(eventlabelspos(ev), yl(2)*0.95, eventlabels{ev}, 'FontSize', 6, 'HorizontalAlignment', eventlabelsalign{ev}); end
        
        
        % ********************
        % Correlations
        
        % Data
        coeffs = paselected.(limbs{b}).(analysis).(toks{1}{2}).coeff(:, str2double(toks{1}{3}));
        normpc2 = paquantl.(limbs{b}).(analysis).(toks{1}{2}).normpc2(:, str2double(toks{1}{3}));
        
        % Plot
        subplot(2, 1, 2);
        yyaxis left;
        plot(0:100, coeffs,'k');
        ylabel('PC coefficient');
        set(gca, 'YColor', 'k');
        yyaxis right;
        hold on;
        plot(0:100, normpc2, 'k--');
        for ev=2:5, xline(events(ev), ':'); end
        hold off;
        ylabel({'Explained variance', ''});
        set(gca, 'YColor', 'k');        
        xlim([0 100]);
        xlabel('% of step-down-and-pivot task');
        %legend('PC coeff', 'Explained var', '', '', '', '');
        

        % ********************
        % Export
        
        % Output folder
        figdir = fullfile(outpath, 'figures');
        if ~exist(figdir, 'dir'), mkdir(figdir); end
           
        % save figures and close
        saveas(gcf, fullfile(figdir, [limbs{b} '_' final.(limbs{b}).labels{f} '.fig']));
        saveas(gcf, fullfile(figdir, [limbs{b} '_' final.(limbs{b}).labels{f} '.png']));
        saveas(gcf, fullfile(figdir, [limbs{b} '_' final.(limbs{b}).labels{f} '.pdf']));
        close(gcf);
        
    end
    

fprintf('------------------------------------------------\n');

end
   



