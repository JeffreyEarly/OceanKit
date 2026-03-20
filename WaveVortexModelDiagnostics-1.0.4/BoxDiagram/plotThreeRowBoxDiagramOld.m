function hFig = plotThreeRowBoxDiagram(row1,row2,row3, arrows, options)
%PLOTBOXDIAGRAM  Lay out three rows of Box objects and scale arrow widths
%
%   h = plotBoxDiagram(row1,row2,row3,arrows,Name,Value,...)  positions the
%   boxes, infers appropriate LineWidth values from each Arrow.Magnitude,
%   then draws everything.

    arguments
        row1 Box
        row2 Box
        row3 Box
        arrows Arrow
        options.RowLabels (1,3) string = ["Sources" "Reservoirs" "Sinks"]
        options.RowSublabels (1,3) string = strings(1,3)
        options.BoxSize (1,2) double {mustBePositive} = [2 1]
        options.Dx double {mustBeNonnegative} = 1
        options.Dy double {mustBeNonnegative} = 2
        options.Origin (1,2) double = [0 9]
        options.Title string = ""
        options.BackgroundColor (1,3) double = [1 1 1]
        options.BaseLineWidth double {mustBePositive} = 1.5
        options.visible = "on"
    end

    rows   = {row1,row2,row3}; nRows = 3;
    w      = options.BoxSize(1); h = options.BoxSize(2);
    dx     = options.Dx;        dy = options.Dy;
    x0     = options.Origin(1); y0 = options.Origin(2);

    maxN   = max(cellfun(@numel,rows));

    layout_style = "fixed-spacing-centered";
    layout_style = "max-spacing";
    if layout_style == "fixed-spacing-centered"
        % Assign positions to boxes
        for r = 1:nRows
            thisRow = rows{r}; n = numel(thisRow);
            xs = x0 + ((maxN - n)*(w + dx))/2 + (0:n-1)*(w + dx);
            y  = y0 - (r-1)*(h + dy);
            for k = 1:n
                thisRow(k).Size = [w h];
                thisRow(k).setPosition([xs(k) y]);
            end
        end
    elseif layout_style == "max-spacing"
        maxWidth = (maxN-1)*(w + dx) + w;
        for r = 1:nRows
            thisRow = rows{r}; n = numel(thisRow);
            if n > 1
                dx = (maxWidth - w)/(n-1)-w;
                xs = x0 + (0:n-1)*(w + dx);
            else
                xs = x0 + (maxWidth-w)/2;
            end    
            y  = y0 - (r-1)*(h + dy);
            for k = 1:n
                thisRow(k).Size = [w h];
                thisRow(k).setPosition([xs(k) y]);
            end
        end
    end



    % Determine line‑width scale from Magnitude field
    mags   = [arrows.Magnitude];
    maxMag = max(mags);
    for a = arrows
        a.LineWidth = options.BaseLineWidth + 4*(a.Magnitude/maxMag);
    end

    % Figure & axes
    hFig = figure('Color',options.BackgroundColor,'Position',[100 100 800 600],Visible=options.visible);
    ax   = axes('Parent',hFig); hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');

    % Draw boxes
    for r = 1:nRows, arrayfun(@(b) b.draw(ax), rows{r}); end

    % Draw arrows
    arrayfun(@(a) a.draw(ax), arrows);

    % Row labels
    lblOffset = 0.8;
    for r = 1:nRows
        yMid  = rows{r}(1).Position(2) + h/2;
        xLeft = x0 - lblOffset;
        text(ax,xLeft,yMid,options.RowLabels(r),'FontSize',16,'FontWeight','bold', ...
            'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
        if all(strlength(options.RowSublabels))
            text(ax,xLeft + 0.3,yMid,options.RowSublabels(r),'FontSize',12,'FontWeight','normal', ...
                'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle');
        end
    end

    if strlength(options.Title)>0
        title(ax,options.Title,'FontWeight','bold', 'Interpreter', 'none');
    end
end