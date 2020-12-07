function ImBat_ProbSuffixTree(out_markov,ord2use);

%compute probabalisit suffix tree for output to cytoscape

%Uses Jeff's code..

% WAL3
%11/26/2020


letters = char(out_markov.VA + 64);
clear DATA;
DATA{1} = letters';
% DATA{1}= convertCharsToStrings(letters);



[F_MAT ALPHABET N PI]=pst_build_trans_mat(DATA,ord2use)

TREE=pst_learn(F_MAT,ALPHABET,N,ord2use);
%pst_export_to_cytoscape(TREE,ALPHABET)

syls = ALPHABET; %convertCharsToStrings(ALPHABET)';
arrange_tree
%[hndls, frqs, degrees] = create_fig7_pies(syls,ALPHABET,TREE);

