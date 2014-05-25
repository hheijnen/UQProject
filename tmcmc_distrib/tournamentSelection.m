%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: tournamentSelection.m
%
% Description: Perform tournament selection without replacement and return
% the matingpool (index of winners of the tournament)
%
% @param tournSize is the tournament size (note that population size should
% be a multiplicant of the tournament size)
%
% @return winners is an array with the index of the individuals that won
% the tournaments.
%
% Author: Kumara Sastry
%
% Date: March 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function winners = tournamentSelection(poolsize,tournSize, fitness)

n = length(fitness);
winners = [];



% Repeat the process "tournSize" times
for i = 1:1:poolsize
    
    %Create a random set of competitors
    shuffleOrder = randperm(n,tournSize);
    competitors = shuffleOrder';
    
    %The winner is the competitor with best fitness
    %[winFit, winID] = max(fitness(competitors),[],2);
    
    winID=find(mnrnd(1,fitness(competitors),1));
    
   % idMap = (0:tournSize-1)*n/tournSize;
    %idMap1 = idMap(winID) + (1:size(competitors,1));
    winners = [winners; competitors(winID)];
end