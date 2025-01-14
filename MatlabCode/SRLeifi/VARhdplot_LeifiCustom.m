function VARhdplot_LeifiCustom(HD,VARopt)
% =======================================================================
% Plot the HD shocks computed with VARhd
% =======================================================================
% VARhdplot(HD,VARopt)
% -----------------------------------------------------------------------
% INPUT
%   - HD: structure from VARhd
%   - VARopt: options of the VAR (see VARopt from VARmodel)
% -----------------------------------------------------------------------
% EXAMPLE
%   - See VARToolbox_Code.m in "../Primer/"
% =======================================================================
% VAR Toolbox 3.0
% Ambrogio Cesa-Bianchi
% ambrogiocesabianchi@gmail.com
% March 2012. Updated April 2021
% -----------------------------------------------------------------------


%% Check inputs
%==========================================================================
if ~exist('VARopt','var')
    error('You need to provide VAR options (VARopt from VARmodel)');
end
% If there is VARopt check that vnames is not empty
vnames = VARopt.vnames;
if isempty(vnames)
    error('You need to add label for endogenous variables in VARopt');
end
% Define shock names
if isempty(VARopt.snames)
    snames = VARopt.vnames;
else
    snames = VARopt.snames;
end

%% Check inputs and define some parameters
%==========================================================================
filename = [VARopt.figname 'HD_'];
quality = VARopt.quality;
suptitle = VARopt.suptitle;
pick = VARopt.pick;

% Initialize HD matrix
[nsteps, nvars, nshocks] = size(HD.shock);

% If one shock is chosen, set the right value for nshocks
if pick<0 || pick>nvars
    error('The selected shock is non valid')
else
    if pick==0
        pick=1;
    else
        nshocks = pick;
    end
end


%% Plot
%==========================================================================
%FigSize(VARopt.FigSize(1),VARopt.FigSize(2))
for ii=pick:nvars
    %H = AreaPlot(squeeze(HD.shock(:,:,ii))); hold on; 
    H = BarPlot(squeeze(HD.shock(:,:,ii))); hold on; 
    h = plot(sum(squeeze(HD.shock(:,:,ii)),2),'-k','LineWidth',2);
    
    if ~isempty(VARopt.firstdate); DatesPlot(VARopt.firstdate,nsteps,8,VARopt.frequency); end
    % Það þarf eitthvað að skoða þessar dagsetingar og stilla 
    % xticks([2:2:136])
    %xticks([2:(4*2):136]);lab = [1990:(1*2):2023.5]; xticklabels(lab)

    xlim([1 nsteps]);
	set(gca,'Layer','top');
    title([vnames{ii}], 'FontWeight','bold','FontSize',10); 
    % Save
    FigName = [filename num2str(ii)];
    if quality 
        if suptitle==1
            Alphabet = char('a'+(1:nvars)-1);
%            SupTitle([Alphabet(jj) ') HD of '  vnames{ii}])
            tempSubtitle = convertCharsToStrings( [Alphabet(ii) ') HD of '  vnames{ii}]);
            SupTitle(tempSubtitle)
        end
        opt = LegOption; opt.handle = [H(1,:) h];
        LegSubplot([snames {'Data'}],opt);
        set(gcf, 'Color', 'w');
        
        export_fig(FigName,'-pdf','-painters');

        %cd ExportFig/;
        %FigNameWithPath = ['../Figures/' FigName];
        %export_fig(FigNameWithPath,'-pdf','-painters')
        %cd ..
    else
        legend([H(1,:) h],[vsnames {'Data'}])
        print('-dpdf','-r100',FigName);
    end
    clf('reset');
end

close all