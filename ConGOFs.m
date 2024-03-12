%Ethel data
units = [1,1.668100537200059,2.782559402207125,4.641588833612778,7.742636826811269,12.915496650148840,21.544346900318835,35.938136638046274,59.948425031894090,100];  %defines X axis units for contrast tuning curves
%File #s. If no LED file, put '0' for YLED and only use NLED.
NLED=2;
YLED=2;

weights = ones(1,10); %Not sure what this is-AJM. %Same :( -SM
a=1;
    i=12;
    fitRates=data(NLED).sumContrast(:,i)';
  fitSterr = ones(1,10);
    [params{a},error{a}] = fit_nakarushtsat(units,fitRates,weights,1);  % modified naka rushton fit to accomodate supersaturating CRFs (Peirce 2007)
    ConFits{a} = nakarushtsat(units,params{a});
    [Powfit,gof(1)] = fit(units',fitRates','power2');
    figure
   errorbar(units, fitRates, fitSterr, '.k');
    hold on
    plot(Powfit,'k');
    xlabel('Contrast (%)');
    ylabel('Spikes/s');
    xlim([0 100]);
    set(gca,'XScale','log');
    title(['Contrast fit for unit ', num2str(i)]);
    hold on
 i=19;
        fitRatesYLED=data(YLED).sumContrast(:,i)';
        fitSterrYLED=ones(1,10);
        [paramsYLED{a},errorYLED{a}] = fit_nakarushtsat(units,fitRatesYLED,weights,1);  % modified naka rushton fit to accomodate supersaturating CRFs (Peirce 2007)
        ConFitsYLED{a} = nakarushtsat(units,paramsYLED{a});
        [PowfitYLED,gof(2)] = fit(units',fitRatesYLED','power2');
        errorbar(units, fitRatesYLED, fitSterrYLED, '.b');
        hold on
        plot(PowfitYLED,'b');

% %Davis data
% units = [0.5,1,1.668100537200059,2.782559402207125,4.641588833612778,7.742636826811269,12.915496650148840,21.544346900318835,35.938136638046274,59.948425031894090,100];  %defines X axis units for contrast tuning curves
% %File #s. If no LED file, put '0' for YLED and only use NLED.
% NLED=1;
% YLED=1;
% 
% weights = ones(1,11); %Not sure what this is-AJM. %Same :( -SM
% a=1;
% 
%     fitRates=data(NLED).con.ave1;
%   fitSterr = ones(1,11)*data(NLED).con.spont1;
%     [params{a},error{a}] = fit_nakarushtsat(units,fitRates,weights,1);  % modified naka rushton fit to accomodate supersaturating CRFs (Peirce 2007)
%     ConFits{a} = nakarushtsat(units,params{a});
%     [Powfit,gof(1)] = fit(units',fitRates','power2');
%     figure
%    errorbar(units, fitRates, fitSterr, '.k');
%     hold on
%     plot(Powfit,'k');
%     xlabel('Contrast (%)');
%     ylabel('Spikes/s');
%     xlim([0 100]);
%     set(gca,'XScale','log');
%     title(['Contrast fit for unit ', num2str(i)]);
%     hold on
% 
%         fitRatesYLED=data(YLED).con.ave3;
%         fitSterrYLED=ones(1,11)*data(NLED).con.spont1;
%         [paramsYLED{a},errorYLED{a}] = fit_nakarushtsat(units,fitRatesYLED,weights,1);  % modified naka rushton fit to accomodate supersaturating CRFs (Peirce 2007)
%         ConFitsYLED{a} = nakarushtsat(units,paramsYLED{a});
%         [PowfitYLED,gof(2)] = fit(units',fitRatesYLED','power2');
%         errorbar(units, fitRatesYLED, fitSterrYLED, '.b');
%         hold on
%         plot(PowfitYLED,'b');
