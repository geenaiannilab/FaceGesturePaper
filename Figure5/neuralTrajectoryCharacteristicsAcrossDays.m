% written GI updated 220104
%
% PUPRPOSE: Calculate reduced neural "state space"
%   where each point in reduced neural state-space represents a BIN in TIME

% velocity dimensions are timepoint x cumulative dimensionality x iteration

clear all;
set(0,'defaultAxesFontSize',24)
set(0,'defaultAxesFontWeight','bold')

subject2analyze ={'combined'};
popType = 'unbalanced';
saveFlag = false; 
workdir = (['/Users/geena/Dropbox/PhD/SUAinfo/Pseudopopulations/' popType '/']);
bhvs = {'Thr','LS','Chew'};
dim2test = [12];

% load data
data = load([workdir 'neuralTrajectory_' subject2analyze{:} '.mat']);

% extract parameters
nIterations = data.nIterations;
arrayList = fieldnames(data.velChewTotal);
win = data.win;
tmin = data.tmin2take;
tmax = data.tmax2take;
taxis = -tmin:win:tmax;  

colorArray = [0.4940 0.1840 0.5560;0.6350 0.0780 0.1840;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250];

%%%%%%%%%%%%%%%%%%
for ss = 1:length(subject2analyze)

    thisSubj = subject2analyze{ss};

    for aa = 1:length(arrayList)

        distanceThrVLS.(arrayList{aa}).mean = mean(data.distanceThrVLS.(arrayList{aa}).data,3);
        distanceThrVLS.(arrayList{aa}).sem = (std(data.distanceThrVLS.(arrayList{aa}).data,[],3)) ./ sqrt(nIterations);
 
        distanceLSVCh.(arrayList{aa}).mean = mean(data.distanceLSVCh.(arrayList{aa}).data,3);
        distanceLSVCh.(arrayList{aa}).sem = (std(data.distanceLSVCh.(arrayList{aa}).data,[],3)) ./ sqrt(nIterations);

        distanceThrVCh.(arrayList{aa}).mean = mean(data.distanceThrVCh.(arrayList{aa}).data,3);
        distanceThrVCh.(arrayList{aa}).sem = (std(data.distanceThrVCh.(arrayList{aa}).data,[],3)) ./ sqrt(nIterations);

        thrVelocity.(arrayList{aa}).mean = mean(data.velThrTotal.(arrayList{aa}).data,3);
        thrVelocity.(arrayList{aa}).sem = (std(data.velThrTotal.(arrayList{aa}).data,[],3)) ./ sqrt(nIterations);

        lsVelocity.(arrayList{aa}).mean = mean(data.velLSTotal.(arrayList{aa}).data,3);
        lsVelocity.(arrayList{aa}).sem = (std(data.velLSTotal.(arrayList{aa}).data,[],3)) ./ sqrt(nIterations);

        chewVelocity.(arrayList{aa}).mean = mean(data.velChewTotal.(arrayList{aa}).data,3);
        chewVelocity.(arrayList{aa}).sem = (std(data.velChewTotal.(arrayList{aa}).data,[],3)) ./ sqrt(nIterations);

        %allBhvVelocity(arrayList{aa}).mean = mean(thrVelocity.(arrayList{aa}).mean, lsVelocity.(arrayList{aa}).mean, chewVelocity.(arrayList{aa}).mean );
    end

    fig1 = figure('Position', [476    83   975   783]);
    for aa = 1:length(arrayList)
        subplot(2,2,aa)
        errorbar(taxis(1:end-1), thrVelocity.(arrayList{aa}).mean(:,end), thrVelocity.(arrayList{aa}).sem(:,end) ,'Color','r','linew',2); hold on
        errorbar(taxis(1:end-1), lsVelocity.(arrayList{aa}).mean(:,end), lsVelocity.(arrayList{aa}).sem(:,end),'Color','b','linew',2);
        errorbar(taxis(1:end-1), chewVelocity.(arrayList{aa}).mean(:,end), chewVelocity.(arrayList{aa}).sem(:,end) ,'Color','g','linew',2); 
        avgTraj = mean([thrVelocity.(arrayList{aa}).mean(:,end), lsVelocity.(arrayList{aa}).mean(:,end), chewVelocity.(arrayList{aa}).mean(:,end)],2);

        plot(taxis(1:end-1), avgTraj,'--k','linew',5);
        ylabel('Speed, au'); xlabel('Time, s'); ylim([0 3]); xlim([-0.8 0.8]);
        title([arrayList{aa}])
    end
    sgtitle('Neural Trajectory Velocities, All Dimensions','FontSize',28,'FontWeight','bold');
    if saveFlag
        saveas(fig1, '/Users/geena/Dropbox/PhD/Manuscript/figures/neuralTraj_byExp_allDims.fig');
        saveas(fig1, '/Users/geena/Dropbox/PhD/Manuscript/figures/neuralTraj_byExp_allDims.png');
    end

    fig2 = figure('Position', [476    83   975   783]);
    for aa = length(arrayList):-1:1
        avgSEM = mean([thrVelocity.(arrayList{aa}).sem(:,end), lsVelocity.(arrayList{aa}).sem(:,end), chewVelocity.(arrayList{aa}).sem(:,end)],2);
        avgTraj2 = mean([thrVelocity.(arrayList{aa}).mean(:,end), lsVelocity.(arrayList{aa}).mean(:,end), chewVelocity.(arrayList{aa}).mean(:,end)],2);
        errorbar(taxis(1:end-1), avgTraj2,avgSEM, 'Color',colorArray(aa,:),'linew',5); hold on 
        ylabel('Speed, au'); xlabel('Time, s'); ylim([0 2]); xlim([-0.8 0.8]); %axes = gca; axes.YTick= [0    0.2000    0.4000    0.6000    0.8000    1.0000];
        l = legend('M3','PMv','M1','S1');
    end
    hold off
    sgtitle('Neural Trajectory Velocities, All Dimensions','FontSize',28,'FontWeight','bold');
    if saveFlag
        saveas(fig2, '/Users/geena/Dropbox/PhD/Manuscript/figures/neuralTraj_Avg_allDims.fig');
        saveas(fig2, '/Users/geena/Dropbox/PhD/Manuscript/figures/neuralTraj_Avg_allDims.png');
    end
%%
    figure; 
    for aa = 1:length(arrayList)
        avgDistance = mean([distanceThrVLS.(arrayList{aa}).mean(:,end) distanceThrVCh.(arrayList{aa}).mean(:,end) distanceLSVCh.(arrayList{aa}).mean(:,end)],2);
        subplot(2,2,aa)
        errorbar(taxis, distanceThrVLS.(arrayList{aa}).mean(:,end), distanceThrVLS.(arrayList{aa}).sem(:,end) ,'Color',[0.3010 0.7450 0.9330],'linew',2); hold on
        errorbar(taxis, distanceThrVCh.(arrayList{aa}).mean(:,end), distanceThrVCh.(arrayList{aa}).sem(:,end) ,'Color',[0.4660 0.6740 0.1880],'linew',2);
        errorbar(taxis, distanceLSVCh.(arrayList{aa}).mean(:,end), distanceLSVCh.(arrayList{aa}).sem(:,end) ,'Color',[0.6350 0.0780 0.1840]	,'linew',2);
        plot(taxis, avgDistance,'k','linew',2); hold off
        ylabel('Distance, au'); xlabel('Time, s'); ylim([5 28 ])
        title([arrayList{aa} '; All Dims'])
        legend('Thr-LS','Thr-Ch','LS-Ch');
    end
    sgtitle([subject2analyze{:} ' All Neural Distances'],'FontSize',20);

    figure;
    for aa = 1:length(arrayList)
        avgDistance = mean([distanceThrVLS.(arrayList{aa}).mean(:,end) distanceThrVCh.(arrayList{aa}).mean(:,end) distanceLSVCh.(arrayList{aa}).mean(:,end)],2);
        avgSEM = sum([distanceThrVLS.(arrayList{aa}).sem(:,end) distanceThrVCh.(arrayList{aa}).sem(:,end) distanceLSVCh.(arrayList{aa}).sem(:,end)],2);
        errorbar(taxis,avgDistance, avgSEM,'Color', colorArray(aa,:),'linew',3.5); hold on
        ylim([0 25]); ylabel('Distance, au'); xlabel('Time, s')

        legend('S1','M1','PMv','M3');
    end
    sgtitle(['Average Intra-Trajectory Distances; All Dimensions'],'FontSize',28);


    %%
    dim2test = [6 8 12 20];
    for dd = 1:length(dim2test)

        figure;
        for aa = 1:length(arrayList)
            avgDistance = mean([distanceThrVLS.(arrayList{aa}).mean(:,dd) distanceThrVCh.(arrayList{aa}).mean(:,dd) distanceLSVCh.(arrayList{aa}).mean(:,dd)],2);
            avgSEM = sum([distanceThrVLS.(arrayList{aa}).sem(:,dd) distanceThrVCh.(arrayList{aa}).sem(:,dd) distanceLSVCh.(arrayList{aa}).sem(:,dd)],2);
            errorbar(taxis,avgDistance, avgSEM,'Color', colorArray(aa,:),'linew',3.5); hold on
            ylim([0 25]); ylabel('Distance, au'); xlabel('Time, s')
            
            legend('S1','M1','PMv','M3');
        end
       sgtitle(['Average Intra-Trajectory Distances; Dims = ' num2str(dim2test(dd)) ],'FontSize',28);

        figure;
        for aa = 1:length(arrayList)
            avgDistance = mean([distanceThrVLS.(arrayList{aa}).mean(:,dd) distanceThrVCh.(arrayList{aa}).mean(:,dd) distanceLSVCh.(arrayList{aa}).mean(:,dd)],2);
            subplot(2,2,aa)
            errorbar(taxis, distanceThrVLS.(arrayList{aa}).mean(:,dd), distanceThrVLS.(arrayList{aa}).sem(:,dd) ,'Color',[0.3010 0.7450 0.9330],'linew',2); hold on
            errorbar(taxis, distanceThrVCh.(arrayList{aa}).mean(:,dd), distanceThrVCh.(arrayList{aa}).sem(:,dd) ,'Color',[0.4660 0.6740 0.1880],'linew',2);
            errorbar(taxis, distanceLSVCh.(arrayList{aa}).mean(:,dd), distanceLSVCh.(arrayList{aa}).sem(:,dd) ,'Color',[0.6350 0.0780 0.1840]	,'linew',2); 
            plot(taxis, avgDistance,'k','linew',2); hold off
            ylim([0 25]); ylabel('Distance, au'); xlabel('Time, s')
            title([arrayList{aa} '; Dims = ' num2str(dim2test(dd))])
            legend('Thr-LS','Thr-Ch','LS-Ch');
        end
        sgtitle([subject2analyze{:} ' All Neural Distances'],'FontSize',20);

        fig3 = figure('Position', [476    83   975   783]);
        for aa = length(arrayList):-1:1
            avgSEM = mean([thrVelocity.(arrayList{aa}).sem(:,dim2test(dd)), lsVelocity.(arrayList{aa}).sem(:,dim2test(dd)), chewVelocity.(arrayList{aa}).sem(:,dim2test(dd))],2);
            avgTraj2 = mean([thrVelocity.(arrayList{aa}).mean(:,dim2test(dd)), lsVelocity.(arrayList{aa}).mean(:,dim2test(dd)), chewVelocity.(arrayList{aa}).mean(:,dim2test(dd))],2);
            errorbar(taxis(1:end-1), avgTraj2,avgSEM, 'Color',colorArray(aa,:),'linew',5); hold on
            ylabel('Speed, au'); xlabel('Time, s'); ylim([0 2]); xlim([-0.8 0.8]); 
            l = legend('M3','PMv','M1','S1');
        end
        hold off
       sgtitle(['Neural Trajectory Velocities, Dim = ' num2str(dim2test(dd))],'FontSize',28,'FontWeight','bold');
       if saveFlag
           saveas(fig3, ['/Users/geena/Dropbox/PhD/Manuscript/figures/neuralTraj_Avg_' num2str(dim2test(dd)) 'D.fig']);
           saveas(fig3, ['/Users/geena/Dropbox/PhD/Manuscript/figures/neuralTraj_Avg_' num2str(dim2test(dd)) 'D.png']);
       end

        fig4 = figure('Position', [476    83   975   783]);
        for aa = 1:length(arrayList)
            subplot(2,2,aa)
            
            errorbar(taxis(1:end-1), thrVelocity.(arrayList{aa}).mean(:,dim2test(dd)), thrVelocity.(arrayList{aa}).sem(:,dim2test(dd)) ,'Color','r','linew',2); hold on
            errorbar(taxis(1:end-1), lsVelocity.(arrayList{aa}).mean(:,dim2test(dd)), lsVelocity.(arrayList{aa}).sem(:,dim2test(dd)),'Color','b','linew',2); 
            errorbar(taxis(1:end-1), chewVelocity.(arrayList{aa}).mean(:,dim2test(dd)), chewVelocity.(arrayList{aa}).sem(:,dim2test(dd)) ,'Color','g','linew',2); 
            avgTraj3 = mean([thrVelocity.(arrayList{aa}).mean(:,dim2test(dd)), lsVelocity.(arrayList{aa}).mean(:,dim2test(dd)), chewVelocity.(arrayList{aa}).mean(:,dim2test(dd))],2);
            plot(taxis(1:end-1), avgTraj3,'--k','linew',5);       
            ylabel('Speed, au'); xlabel('Time, s'); ylim([0 2.2]); xlim([-0.8 0.8]);
            title([arrayList{aa}])

        end
       sgtitle(['Neural Trajectory Speeds, ' num2str(dim2test(dd)) '-D'],'FontSize',28,'FontWeight','bold');
       if saveFlag
           saveas(fig4, ['/Users/geena/Dropbox/PhD/Manuscript/figures/neuralTraj_byExp_' num2str(dim2test(dd)) 'D.fig']);
           saveas(fig4, ['/Users/geena/Dropbox/PhD/Manuscript/figures/neuralTraj_byExp_' num2str(dim2test(dd)) 'D.png']);
       end
    end % end dimensions to test


    for aa = 1:length(arrayList)

        fig5 = figure('Position', [1    85   384   781]);
        dim2test = [8 12 20];
        
        for dd = 1:length(dim2test)
        
            subplot(3,1,dd)
            errorbar(taxis(1:end-1), thrVelocity.(arrayList{aa}).mean(:,dd), thrVelocity.(arrayList{aa}).sem(:,dd) ,'Color','r','linew',2); hold on
            errorbar(taxis(1:end-1), lsVelocity.(arrayList{aa}).mean(:,dd), lsVelocity.(arrayList{aa}).sem(:,dd),'Color','b','linew',2); 
            errorbar(taxis(1:end-1), chewVelocity.(arrayList{aa}).mean(:,dd), chewVelocity.(arrayList{aa}).sem(:,dd) ,'Color','g','linew',2); hold off
            ylim([0 2.2]); xlim([-0.8 0.8]); ylabel('Speed, au'); xlabel('Time, s')
            title([arrayList{aa} '; Dims = ' num2str(dim2test(dd))])
    
        end
        if saveFlag
            saveas(fig5, ['/Users/geena/Dropbox/PhD/Manuscript/figures/neuralTrajVelocities_' arrayList{aa} '.fig']);
            saveas(fig5, ['/Users/geena/Dropbox/PhD/Manuscript/figures/neuralTrajVelocities_' arrayList{aa} '.png']);
        end
    end
%%



end