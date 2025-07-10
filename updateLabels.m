function updatedLabels = updateLabels(labels, uniqueLabels, newLabels)

    uniqueLabels = uniqueLabels(:);
    newLabels = newLabels(:);
    
    [~, loc] = ismember(labels, uniqueLabels);
        
    if any(loc == 0)
        error('Some labels in ''labels'' do not exist in ''uniqueLabels''.');
    end
        
    updatedLabels = newLabels(loc);
end