echo "Run Target"

#python make_target_peps_general.py /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Mutator_pep/fasta /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/reference/AaegL3.fa target 
#qsub target_aedes_MULE_tpase.sh
#qsub target_eukaryote_MULE_DDE_domain.sh

echo "Generate nonredundant BEDfile and sequence files"
#python make_BEDfile_all_protein_hits.py /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/all_pep-hits.bed /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/reference/AaegL3.fa /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/nonredundant_pep_union-10kb.fa /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/nonredundant_pep_union-hits.fa

echo "Seperate flank sequences to individual files"
#python split_fasta.py /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/nonredundant_pep_union-10kb.fa /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/split

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

echo "Generate sequence files of each clusterwith >= 2 seqs and 1 only"
perl /opt/cd-hit/4.6.1/bin/make_multi_seq.pl /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/nonredundant_pep_union-hits.fa /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/nonredundant_pep_union-hits_c90.clstr /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/Autonomous/Target/nonredundant_pep_union-hits_c90_seq 2

perl make_multi_seq.pl Target/nonredundant_pep_union-hits.fa Target/nonredundant_pep_union-hits_c90.clstr Target/nonredundant_pep_union-hits_c90_multi 2
perl make_multi_seq.pl Target/nonredundant_pep_union-hits.fa Target/nonredundant_pep_union-hits_c90.clstr Target/nonredundant_pep_union-hits_c90_single 1

echo "Rename protein hit sequences to include predicted pee match results"
python rename_hit-flank-files_by_peps-to-repbase_results.py Target/predicted_to_repbase_banshee_full.txt Target/predicted_to_repbase_banshee_retros.txt Target/predicted_to_repbase_banshee_heli-mav.txt Target/predicted_to_repbase_banshee_no-matches.txt Target/nonredundant_pep_union-hits.fa Target/nonredundant_pep_union-10kb.fa

echo "move clusters that have mojority predictions being retro, helitron or maverick. Also make new sequence files cleaned of the unwanted sequences"


echo "summary family predicted in fasta"
python summarize_protein_matches.py Target/nonredundant_pep_union-hits_match-info.fa Target/nonredundant_pep_union-hits_match-info.fa.summary

