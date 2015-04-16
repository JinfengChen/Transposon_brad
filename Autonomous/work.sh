##Start annotation
echo "Run Target"
python make_target_peps_general.py /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Mutator_pep/fasta /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/reference/AaegL3.fa target 
qsub target_aedes_MULE_tpase.sh
qsub target_eukaryote_MULE_DDE_domain.sh

echo "Generate nonredundant BEDfile and sequence files"
python make_BEDfile_all_protein_hits.py /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/all_pep-hits.bed /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/reference/AaegL3.fa /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/nonredundant_pep_union-10kb.fa /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/nonredundant_pep_union-hits.fa

echo "Seperate flank sequences to individual files"
python split_fasta.py /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/nonredundant_pep_union-10kb.fa /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/split

echo "Generate Augustus prediction"
python make_augustus_directory.py /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/split aedes fly
qsub augustus_all.sh
cat /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/split/*.gff > /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/augustus.gff

echo "Parse prediction from those in initial hit region"
python parse_augustus_for_proteins_in_range.py /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/augustus.gff

echo "Generate genscan prediction"
python make_genscan.py /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/split *.fa human
qsub genscan_all.sh
cat /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/split/*.out > /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/genscan.out

echo "parse genscan predictions"
python parse_genscan.py /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/genscan.out

echo "Concatenate parsed augustus and genscan prediction proteins"
cat /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/augustus_protein.fa /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/genscan_protein.fa > /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/nrprotein_augustus_genscan.fa

echo "Search wanted prediction against Repbase proteins with glsearch. Parse results"
python make_ggsearch_onefile.py /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/nrprotein_augustus_genscan.fa /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/repbase/repbase18.08-aa_bansee.fa
qsub ggsearch-nrprotein_augustus_genscan-to_repbase18.08-aa_bansee.sh
module load perl/5.10.1
perl match_predicted-peps_glsearch.pl /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/nrprotein_augustus_genscan_ggsearch-to_repbase18.08-aa_bansee.out predicted_to_repbase_banshee
 
echo "Cluster nonredundant hit sequence"
python make_cdhit-prot_sh.py /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/nonredundant_pep_union-hits.fa
qsub cd-hit_nonredundant_pep_union-hits-90.sh
python make_multi_seq_sum.py Target/nonredundant_pep_union-hits_c90.clstr
python make_multi_seq_sum.py Target/nonredundant_pep_union-hits_c80.clstr

echo "Generate sequence files of each clusterwith >= 2 seqs and 1 only"
perl make_multi_seq.pl Target/nonredundant_pep_union-10kb.fa Target/nonredundant_pep_union-hits_c90.clstr Target/nonredundant_pep_union-10kb_c90_multi 2
perl make_multi_seq.pl Target/nonredundant_pep_union-10kb.fa Target/nonredundant_pep_union-hits_c90.clstr Target/nonredundant_pep_union-10kb_c90_single 1
perl make_multi_seq.pl Target/nonredundant_pep_union-hits.fa Target/nonredundant_pep_union-hits_c90.clstr Target/nonredundant_pep_union-hits_c90_multi 2
perl make_multi_seq.pl Target/nonredundant_pep_union-hits.fa Target/nonredundant_pep_union-hits_c90.clstr Target/nonredundant_pep_union-hits_c90_single 1

perl make_multi_seq.pl Target/nonredundant_pep_union-10kb.fa Target/nonredundant_pep_union-hits_c80.clstr Target/nonredundant_pep_union-10kb_c80_multi 2
perl make_multi_seq.pl Target/nonredundant_pep_union-10kb.fa Target/nonredundant_pep_union-hits_c80.clstr Target/nonredundant_pep_union-10kb_c80_single 1

##Can stop here for MULE project for Kun, just use activeTE to find TSD and TIR

##clean data when doing annotation for all DNA transposon
echo "Rename protein hit sequences to include predicted pee match results"
python rename_hit-flank-files_by_peps-to-repbase_results.py Target/predicted_to_repbase_banshee_full.txt Target/predicted_to_repbase_banshee_retros.txt Target/predicted_to_repbase_banshee_heli-mav.txt Target/predicted_to_repbase_banshee_no-matches.txt Target/nonredundant_pep_union-hits.fa Target/nonredundant_pep_union-10kb.fa

echo "move clusters that have mojority predictions being retro, helitron or maverick. Also make new sequence files cleaned of the unwanted sequences"


echo "summary family predicted in fasta"
python summarize_protein_matches.py Target/nonredundant_pep_union-hits_match-info.fa Target/nonredundant_pep_union-hits_match-info.fa.summary

##activeTE for autonomous cluster
echo "Trim 10 kb flanks down to 5kb (will use 10kb if need, i.e. flanks are too conserved by aTE analysis)"
python trim_flanks_directory.py Target/nonredundant_pep_union-10kb_c90_multi/ 5000
python trim_flanks_directory.py Target/nonredundant_pep_union-10kb_c80_multi/ 5000

echo "Run multiple sequence alignments"
python make_mafft_pep-cluster_split_one_sh.py Target/nonredundant_pep_union-10kb_c90_multi/ '*trimmed-5000.fa' aedes
python make_mafft_pep-cluster_split_one_sh.py Target/nonredundant_pep_union-10kb_c80_multi/ '*trimmed-5000.fa' aedes
qsub -q js mafft_aedes.sh

echo "Analyzing MSA files with activeTE"
python make_activeTE-pep-msa-one.py Target/nonredundant_pep_union-10kb_c90_multi/ '*.msa' c90
python make_activeTE-pep-msa-one.py Target/nonredundant_pep_union-10kb_c80_multi/ '*.msa' c80
qsub aTE-pep_c80.sh

echo "Initial summarizing of aTE"
python list_aTE_good_bad_mixed_pep.py Target/nonredundant_pep_union-10kb_c90_multi/aTE Target/nonredundant_pep_union-10kb_c90_multi/aTE_sum aedes_pep
python list_aTE_good_bad_mixed_pep.py Target/nonredundant_pep_union-10kb_c80_multi/aTE Target/nonredundant_pep_union-10kb_c80_multi/aTE_sum aedes_pep

echo "Recombine results of split groups"
python recombine_aTE-out3.py Target/nonredundant_pep_union-10kb_c90_multi/aTE_sumaedes_pep.good Target/nonredundant_pep_union-10kb_c90_multi/aTE/ Target/nonredundant_pep_union-10kb_c90_multi/aTE_recombine aedes_pep


