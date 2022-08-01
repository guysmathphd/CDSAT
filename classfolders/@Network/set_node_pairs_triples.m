function obj = set_node_pairs_triples(obj)

k = General.load_var(fullfile(obj.path,'degree_vector.mat'));
[~,I] = sort(k,'descend');

doubles = nchoosek(I,2);
disp(['size of doubles is: ' num2str(size(doubles,1)) ' by ' num2str(size(doubles,2))]);
General.save_var(doubles,fullfile(obj.path,'node_combs'),'doubles');

subset_size = 1000;
triples = nchoosek(I(1:subset_size),3);
disp(['size of triples is: ' num2str(size(triples,1)) ' by ' num2str(size(triples,2))]);
General.save_var(triples,fullfile(obj.path,'node_combs'),'triples');

end