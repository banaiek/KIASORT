function applyColorScheme(h, parentBg)
brightness = mean(parentBg);
if brightness < 0.5
    textColor = [1 1 1];
else
    textColor = [0 0 0];
end
if isprop(h,'BackgroundColor')
    h.BackgroundColor = parentBg;
end
if isprop(h,'FontColor')
    h.FontColor = textColor;
end
if isprop(h,'ForegroundColor')
    h.ForegroundColor = textColor;
end
if isprop(h,'Children')
    kids = h.Children;
    for k = 1:numel(kids)
        applyColorScheme(kids(k), parentBg);
    end
end
end
