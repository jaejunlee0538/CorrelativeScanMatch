function l=prob_to_log_odds(p)
% Convert proability values p to the corresponding log odds l.
% p could be a scalar or a matrix.

l = p ./ (1 - p);
l = log(l);
% TODO: compute l.


end
