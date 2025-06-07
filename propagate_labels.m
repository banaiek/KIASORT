function finalVec = propagate_labels(firstVec, finalVec)

    changed   = firstVec ~= finalVec;      % logical index of changes
    from      = firstVec(changed);         % original labels that moved
    to        = finalVec(changed);         % where they moved to

    [~, ia]   = unique(from, 'stable');
    from      = from(ia);
    to        = to(ia);

    while true
        [tf, loc]   = ismember(finalVec, from);   
        newFinal    = finalVec;                  
        newFinal(tf) = to(loc(tf));               

        if isequal(newFinal, finalVec)            % converged
            break
        end
        finalVec = newFinal;                     
    end
end