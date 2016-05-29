function h = normalizeHeading(heading)
    if heading > 0
        h = mod(heading, 2*pi);
        if h > pi
            h = -2*pi + h;
        end
    else
        h = -mod(-heading, 2*pi);
        if h < -pi
            h = 2*pi + h;
        end
    end
end

