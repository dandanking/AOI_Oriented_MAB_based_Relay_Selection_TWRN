 function [] = display_results_individual_performance(aoiEvolution, aoi_benchmark_MMRS, aoi_benchmark_MS,...
        tx_power_SR, tx_power_RD, timesArmHasBeenPlayed )
% display_results_individual_performance displays the results of experiments for
% individual performance of each policy 
%   INPUT: 
%       * aoiEvolution - array containing the average aoi experienced
%        for each of the experiment iterations
%       * timesArmHasBeenPlayed - array containing information about the
%       number of plays that WN has done for each action

    % Load constatns
    load('constants.mat')
    avg_txpower_SR=zeros(2, totalIterations);
    avg_txpower_RD=zeros(2, totalIterations);
    tx_power_SR_mW = zeros(2, totalIterations);
    tx_power_RD_mW = zeros(2, totalIterations);
    for iter=1:totalIterations
        for i=1:2
            tx_power_SR_mW(i,iter) = db2pow(tx_power_SR(i,iter));
        end
    end
    
    for iter=1:totalIterations
        a=sum(tx_power_SR_mW(:,1:iter),2)./iter;
        for i=1:2
            avg_txpower_SR(i,iter)=a(i);
        end
    end
    for iter=1:totalIterations
        for i=1:2
            if tx_power_RD(i,iter) == 1
                tx_power_RD_mW(i,iter) = db2pow(2).*0.5;
            else
                tx_power_RD_mW(i,iter) = db2pow(2);
            end
        end
    end
    for iter=1:totalIterations
        a=sum(tx_power_RD_mW(:,1:iter),2)./iter;
        for i=1:2
            avg_txpower_RD(i,iter)=a(i);
        end
    end
    benchmarkTP=ones(2, totalIterations).*db2pow(2);
    % Set font type
    set(0,'defaultUicontrolFontName','Times New Roman');
    set(0,'defaultUitableFontName','Times New Roman');
    set(0,'defaultAxesFontName','Times New Roman');
    set(0,'defaultTextFontName','Times New Roman');
    set(0,'defaultUipanelFontName','Times New Roman');

    %% average aoi experienced for each iteration
    fig = figure('pos',[450 400 500 350]);
    axis([1 20 30 70]);
    aoi_per_iteration = aoiEvolution(1:totalIterations);
    benchmark_aoi_MMRS = aoi_benchmark_MMRS(1:totalIterations);
    benchmark_aoi_MS = aoi_benchmark_MS(1:totalIterations);
    plot(1:totalIterations, aoi_per_iteration, 1:totalIterations, benchmark_aoi_MMRS,'r', 1:totalIterations, benchmark_aoi_MS,'g');
    hold on
    axis([1 totalIterations 0 20])
    xlabel( ' iteration', 'fontsize', 24)
    ylabel('average aoi', 'fontsize', 24)   
    
    %% Actions probability
    fig = figure('pos',[450 400 500 350]);
    axis([1 20 30 70]); 
    % Print the preferred action per wlan   
    for i = 1:2
        K = size(timesArmHasBeenPlayed, 2);
        subplot(2,1,i);
        bar(1:K, timesArmHasBeenPlayed(i, :)/totalIterations);
        hold on
        axis([0 K+1 0 1])
        xticks(1:K)
        xticklabels(1:K)
        xlabel( 'Action index', 'fontsize', 20)
        ylabel('Action prob.', 'fontsize', 20)   
    end
    
    %% Average transmission power first-hop
    fig = figure('pos',[450 400 500 350]);
    axis([1 20 30 70]); 
    % Print the average transmission power per wlan   
    for i = 1:2
        subplot(2,1,i);
        plot(1:totalIterations, avg_txpower_SR(i,:), 1:totalIterations, benchmarkTP(i,:),'r');
        hold on
        axis([1 totalIterations -2 3])
        xlabel( 'iteration', 'fontsize', 20)
        ylabel('Average TP', 'fontsize', 20)   
    end
     %% Average transmission power second-hop
    fig = figure('pos',[450 400 500 350]);
    axis([1 20 30 70]); 
    % Print the average transmission power per wlan   
    for i = 1:2
        subplot(2,1,i);
        plot(1:totalIterations, avg_txpower_RD(i,:), 1:totalIterations, benchmarkTP(i,:),'r');
        hold on
        axis([1 totalIterations 0 3])
        xlabel( 'iteration', 'fontsize', 20)
        ylabel('Average TP', 'fontsize', 20)   
    end
end