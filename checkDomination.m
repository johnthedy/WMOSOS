function dom_vector = checkDomination(fitness)
Np = size(fitness,1);
dom_vector = zeros(Np,1);
all_perm = nchoosek(1:Np,2);    % Possible permutations
all_perm = [all_perm; [all_perm(:,2) all_perm(:,1)]];

d = dominates(fitness(all_perm(:,1),:),fitness(all_perm(:,2),:));
dominated_particles = unique(all_perm(d==1,2));
dom_vector(dominated_particles) = 1;
end

% Function that returns 1 if x dominates y and 0 otherwise
function d = dominates(x,y)
d = all(x<y,2) & any(x<y,2);
end