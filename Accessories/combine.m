function combination = combine(small_vec, background_vec, location)
    if isempty(location)
        combination = small_vec;
    else
        combination = repmat(background_vec, size(small_vec, 1), 1);
        combination(:, location) = small_vec;
    end
end