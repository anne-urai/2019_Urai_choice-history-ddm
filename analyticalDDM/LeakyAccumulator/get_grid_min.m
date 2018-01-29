function [max_ps,max_Ls] = get_grid_min(gridSets,gridL,grid_p,set_p,nOpt)

grid_p_rep{1}=grid_p.beta; grid_p_rep{2}=grid_p.sigmaA; grid_p_rep{3}=grid_p.theta; grid_p_rep{4}=grid_p.leak; 
fixed_ps = find(strcmp({set_p.beta,set_p.sigmaA,set_p.theta,set_p.leak},'fi')==1);  % finding out fixed parameters (i.e. those that limit viable across-condition combinations)

max_ps={}; max_Ls=zeros(1,nOpt);
test_ps={}; test_Ls={};

if isempty(fixed_ps)  % if no parameters are fixed across conditions
    test_ps{1,1} = gridSets;
    test_Ls{1,1} = gridL(:,1)*gridL(:,2)';
    for x = 1:size(test_Ls{1},1);
        for y = 1:size(test_Ls{1},2);
            [minval,minpos] = min(max_Ls);
            if test_Ls{1}(x,y) > minval
                max_Ls(minpos) = test_Ls{1}(x,y);
                max_ps{minpos}(1,1:length(test_ps{1}(x,:))) = test_ps{1}(x,:);
                max_ps{minpos}(2,1:length(test_ps{1}(y,:))) = test_ps{1}(y,:);
            end
        end
    end
    
elseif length(fixed_ps)==1   % if one parameter is fixed across conditions
    for p1 = 1:length(grid_p_rep{fixed_ps});
        test_ps{p1,1} = gridSets(find(gridSets(:,fixed_ps)==grid_p_rep{fixed_ps}(p1)),:);  % store useable parameter sets given constraints
        test_Ls{p1,1} = gridL(find(gridSets(:,fixed_ps)==grid_p_rep{fixed_ps}(p1)),1)*...  % calculate full across-condition likelihoods for all possible parameter sets within these constraints
                      gridL(find(gridSets(:,fixed_ps)==grid_p_rep{fixed_ps}(p1)),2)';
    end
    for p1 = 1:size(test_ps,1)  % search through likelihoods of all possible parameter sets and store n best fits
        for x = 1:size(test_Ls{p1},1);
            for y = 1:size(test_Ls{p1},2);
                [minval,minpos] = min(max_Ls);
                if test_Ls{p1}(x,y) > minval
                    max_Ls(minpos) = test_Ls{p1}(x,y);
                    max_ps{minpos}(1,1:length(test_ps{p1}(x,:))) = test_ps{p1}(x,:);
                    max_ps{minpos}(2,1:length(test_ps{p1}(y,:))) = test_ps{p1}(y,:);
                end
            end
        end
    end
    
elseif length(fixed_ps)==2    % if two parameters are fixed across conditions
    for p1 = 1:length(grid_p_rep{fixed_ps(1)});
        for p2 = 1:length(grid_p_rep{fixed_ps(2)});
            test_ps{p1,p2} = gridSets(find(gridSets(:,fixed_ps(1))==grid_p_rep{fixed_ps(1)}(p1) & gridSets(:,fixed_ps(2))==grid_p_rep{fixed_ps(2)}(p2)),:);
            test_Ls{p1,p2} = gridL(find(gridSets(:,fixed_ps(1))==grid_p_rep{fixed_ps(1)}(p1) & gridSets(:,fixed_ps(2))==grid_p_rep{fixed_ps(2)}(p2)),1)*...
                gridL(find(gridSets(:,fixed_ps(1))==grid_p_rep{fixed_ps(1)}(p1) & gridSets(:,fixed_ps(2))==grid_p_rep{fixed_ps(2)}(p2)),2)';
        end
    end
    for p1 = 1:size(test_ps,1)
        for p2 = 1:size(test_ps,2)
            for x = 1:size(test_Ls{p1,p2},1);
                for y = 1:size(test_Ls{p1,p2},2);
                    [minval,minpos] = min(max_Ls);
                    if test_Ls{p1,p2}(x,y) > minval
                        max_Ls(minpos) = test_Ls{p1,p2}(x,y);
                        max_ps{minpos}(1,1:length(test_ps{p1,p2}(x,:))) = test_ps{p1,p2}(x,:);
                        max_ps{minpos}(2,1:length(test_ps{p1,p2}(y,:))) = test_ps{p1,p2}(y,:);
                    end
                end
            end
        end
    end
    
elseif length(fixed_ps)==3    % if three parameters are fixed across conditions
    for p1 = 1:length(grid_p_rep{fixed_ps(1)});
        for p2 = 1:length(grid_p_rep{fixed_ps(2)});
            for p3 = 1:length(grid_p_rep{fixed_ps(3)});
                test_ps{p1,p2,p3} = gridSets(find(gridSets(:,fixed_ps(1))==grid_p_rep{fixed_ps(1)}(p1) & gridSets(:,fixed_ps(2))==grid_p_rep{fixed_ps(2)}(p2) & gridSets(:,fixed_ps(3))==grid_p_rep{fixed_ps(3)}(p3)),:);
                test_Ls{p1,p2,p3} = gridL(find(gridSets(:,fixed_ps(1))==grid_p_rep{fixed_ps(1)}(p1) & gridSets(:,fixed_ps(2))==grid_p_rep{fixed_ps(2)}(p2) & gridSets(:,fixed_ps(3))==grid_p_rep{fixed_ps(3)}(p3)),1)*...
                    gridL(find(gridSets(:,fixed_ps(1))==grid_p_rep{fixed_ps(1)}(p1) & gridSets(:,fixed_ps(2))==grid_p_rep{fixed_ps(2)}(p2) & gridSets(:,fixed_ps(3))==grid_p_rep{fixed_ps(3)}(p3)),2)';
            end
        end
    end
    for p1 = 1:size(test_ps,1)
        for p2 = 1:size(test_ps,2)
            for p3 = 1:size(test_ps,3)
                for x = 1:size(test_Ls{p1,p2,p3},1);
                    for y = 1:size(test_Ls{p1,p2,p3},2);
                        [minval,minpos] = min(max_Ls);
                        if test_Ls{p1,p2,p3}(x,y) > minval
                            max_Ls(minpos) = test_Ls{p1,p2,p3}(x,y);
                            max_ps{minpos}(1,1:length(test_ps{p1,p2,p3}(x,:))) = test_ps{p1,p2,p3}(x,:);
                            max_ps{minpos}(2,1:length(test_ps{p1,p2,p3}(y,:))) = test_ps{p1,p2,p3}(y,:);
                        end
                    end
                end
            end
        end
    end
    
end