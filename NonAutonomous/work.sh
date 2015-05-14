echo "Run repeatmoduler and MITEhunter, merge the resulted repeat elements into one file"
#cp ../MITEhunter/AaegL3_Step8.MH.fa ./
#cp ../RepeatModeler/RM_1251.SatApr42258082015/AaegL3.RM.fa ./
#cat AaegL3_Step8.MH.fa AaegL3.RM.fa > AaegL3.MH_RM.fa
#mv AaegL3*.fa initial_data/

#check output for duplicates and TIRs
echo "Do TIRs searchs, parse output and seperate sequence"
#perl check_query_tirs_ggsearch4.pl initial_data/AaegL3.MH_RM.fa

echo "Run ggsearch on good-tirs files"
#ggsearch36 -E 1e-20 -n -T 2 initial_data/AaegL3.MH_RM/AaegL3.MH_RM_good-tirs.fa initial_data/AaegL3.MH_RM/AaegL3.MH_RM_good-tirs.fa > initial_data/AaegL3.MH_RM/AaegL3.MH_RM_good-tirs-to-self.ggsearch_out

echo "Parse ggsearch output: generate unique list of good tirs"
#perl match_mites_ggsearch_both-ori2.pl initial_data/AaegL3.MH_RM/AaegL3.MH_RM_good-tirs-to-self.ggsearch_out AedesL3 1e-30 same

echo "Separate sequence: generate unique good tirs fasta"
#python separate_sequence_by-list.py initial_data/AaegL3.MH_RM/AedesL3_unique.txt initial_data/AaegL3.MH_RM/AaegL3.MH_RM_good-tirs.fa initial_data/AaegL3.MH_RM/AedesL3_unique.fa
 
echo "Run ggsearch for bad-tirs aganist unique file"
#ggsearch36 -E 1e-20 -n -T 2 initial_data/AaegL3.MH_RM/AaegL3.MH_RM_bad-tirs_trimmed.fa initial_data/AaegL3.MH_RM/AedesL3_unique.fa > initial_data/AaegL3.MH_RM/AaegL3.MH_RM_bad-tirs_trimmed-to-unique.ggsearch_out

echo "Parse ggsearch output: "
#perl match_mites_ggsearch_both-ori.pl initial_data/AaegL3.MH_RM/AaegL3.MH_RM_bad-tirs_trimmed-to-unique.ggsearch_out AedesL3_bad-to-unique 1e-30

echo "Seperate sequence: generate unique list of bad tir that do not matched to good unique tirs"
#python separate_sequence_by-list.py initial_data/AaegL3.MH_RM/AedesL3_bad-to-unique_unique.txt initial_data/AaegL3.MH_RM/AaegL3.MH_RM_bad-tirs_trimmed.fa initial_data/AaegL3.MH_RM/AedesL3_bad-to-unique_unique.fa

echo "prepare target input"
#mkdir initial_data/AaegL3.MH_RM/query
#cat initial_data/AaegL3.MH_RM/initial_data/AaegL3.MH_RM/AedesL3_bad-to-unique_unique.fa initial_data/AaegL3.MH_RM/initial_data/AaegL3.MH_RM/AedesL3_unique.fa initial_data/AaegL3.MH_RM/query/AedesL3_unique_merge.fa

echo "Run target: should split sequence by length, <1000 bp use multiprocess to run multijob (-P 1 -C 16) and >1000 to use multiprocess run single job (-P 16 -C 1)"
#python make_target_nonauto_general.py /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/NonAutonomous/initial_data/AaegL3.MH_RM/query /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/NonAutonomous/Target /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/reference/AaegL3.fa Target_Run
#qsub Target_Run_AedesL3_unique.sh
#qsub Target_Run_AedesL3_bad-to-unique_unique.sh

echo "Run activeTE"
mkdir Target/Group_msa
cp Target/Target_Run_2015_04_30_130853/AedesL3_unique_split*/*.msa Target/Group_msa
python make_activeTE_sh.py Target/Group_msa/ '*.msa' activeTE_test
python multi_sub_N_move.py . '*_under.sh' 100 aTE_submitted/
python multi_sub_N_move.py . '*_split.sh' 100 aTE_submitted/
python list_aTE_good_bad_mixed.py Target/Group_msa Target/Group_msa AedesL3_unique_test_aTE_sum
python categorize_all_aTE_results_by_group.py Target/Group_msa Target/ Group_msa_sum

echo "Run activeTE: unique"
cd Target
mkdir Target_Run_unique_group_split_msa
cp Target_Run_unique_finished/AedesL3_unique_split*/*_under.group1_split.msa Target_Run_unique_group_split_msa/
python make_activeTE_sh.py Target/Target_Run_unique_group_split_msa '*.msa' activeTE_Target_Run_unique_group_split
python multi_sub_N_move.py . '*_split.sh' 200 aTE_submitted/
python list_aTE_good_bad_mixed.py Target/Target_Run_unique_group_split_msa/ Target/Target_Run_unique_group_split_msa/ AedesL3_unique_Target_Run_unique_group_split_sum
python categorize_all_aTE_results_by_group.py Target/Target_Run_unique_group_split_msa Target/Target_Run_unique_group_split_msa AedesL3_unique_Target_Run_unique_group_split_sum

mkdir Target_Run_unique_group_msa
cp Target_Run_unique_finished/AedesL3_unique_split*/*_under.msa Target_Run_unique_group_msa/
python make_activeTE_sh.py Target/Target_Run_unique_group_msa/ '*.msa' activeTE_Target_Run_unique_group
python multi_sub_N_move.py . '*_under.sh' 200 aTE_submitted/
python list_aTE_good_bad_mixed.py Target/Target_Run_unique_group_msa/ Target/Target_Run_unique_group_msa/ AedesL3_unique_Target_Run_unique_group_sum
python categorize_all_aTE_results_by_group.py Target/Target_Run_unique_group_msa Target/Target_Run_unique_group_msa AedesL3_unique_Target_Run_unique_group_sum


echo "Redo bad, notsd and conserved, flanking=300bp"
cat Target/Target_Run_unique_group_*/*.bad Target/Target_Run_unique_group_*/*.no_tsd Target/Target_Run_unique_group_*/*.conserved_flanks > Target/Target_Run_unique_group.redo.list
python split_fasta_id.py initial_data/AaegL3.MH_RM/query/AedesL3_unique.fa
python make_redo_fasta.py Target/Target_Run_unique_group.redo.list initial_data/AaegL3.MH_RM/query/AedesL3_unique.fa
python make_target_nonauto_general_redo.py /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/NonAutonomous/initial_data/AaegL3.MH_RM/query/AedesL3_unique.redo/ /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/NonAutonomous/Target /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/reference/AaegL3.fa Target_Run_unique_redo
qsub Target_Run_unique_redo_AedesL3_unique.redo_1k.sh
qsub Target_Run_unique_redo_AedesL3_unique.redo_1kto3k.sh

echo "Run activeTE: bad to unique"
mkdir Target_Run_bad-to-unique_unique_group_split_msa
cp Target_Run_bad-to-unique_unique_finished/AedesL3_bad-to-unique_unique_split*/*_under.group1_split.msa Target_Run_bad-to-unique_unique_group_split_msa/
python make_activeTE_sh.py Target/Target_Run_bad-to-unique_unique_group_split_msa '*.msa' activeTE_Target_Run_bad-to-unique_unique_group_split
python multi_sub_N_move.py . '*_split.sh' 200 aTE_submitted/
python list_aTE_good_bad_mixed.py Target/Target_Run_bad-to-unique_unique_group_split_msa/ Target/Target_Run_bad-to-unique_unique_group_split_msa/ AedesL3_unique_Target_Run_bad-to-unique_unique_group_split_sum

mkdir Target_Run_bad-to-unique_unique_group_msa
cp Target_Run_bad-to-unique_unique_finished/AedesL3_bad-to-unique_unique_split*/*_under.msa Target_Run_bad-to-unique_unique_group_msa/
python make_activeTE_sh.py Target/Target_Run_bad-to-unique_unique_group_msa '*.msa' activeTE_Target_Run_bad-to-unique_unique_group
python multi_sub_N_move.py . '*_under.sh' 200 aTE_submitted/
python list_aTE_good_bad_mixed.py Target/Target_Run_bad-to-unique_unique_group_msa/ Target/Target_Run_bad-to-unique_unique_group_msa/ AedesL3_unique_Target_Run_bad-to-unique_unique_group_sum

echo "Redo bad, notsd and conserved, flanking=300bp"
cat Target/Target_Run_bad-to-unique_unique_group_*/*.bad Target/Target_Run_bad-to-unique_unique_group_*/*.no_tsd Target/Target_Run_bad-to-unique_unique_group_*/*.conserved_flanks > Target/Target_Run_bad-to-unique_group.redo.list
python split_fasta_id.py initial_data/AaegL3.MH_RM/query/AedesL3_bad-to-unique_unique.fa
python make_redo_fasta.py Target/Target_Run_bad-to-unique_group.redo.list initial_data/AaegL3.MH_RM/query/AedesL3_bad-to-unique_unique.fa
python make_target_nonauto_general_redo.py /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/NonAutonomous/initial_data/AaegL3.MH_RM/query/AedesL3_bad-to-unique_unique.redo/ /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/NonAutonomous/Target /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/reference/AaegL3.fa Target_Run_bad-to-unique_redo


echo "TSD and TIR length"
cat Target/*_group_*msa/*.flank_filter-1.2_under*/*.element_info > good.list
cut -f 8 good.list | sort -n
awk '{print length($4)}' good.list | awk '$1>50' | wc -l


echo "MULE related query only"
python make_target_nonauto_general_redo.py /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/NonAutonomous/initial_data_MULE/Merge/query/ /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/NonAutonomous/Target /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/reference/AaegL3.fa Target_Run
python make_target_nonauto_general_redo.py /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/NonAutonomous/initial_data_MULE/Merge/query_noMULE /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/NonAutonomous/Target /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/reference/AaegL3.fa Target_Run
qsub -q highmem Target_Run_AedesL3_MULE.sh #2 cpu and 40g mem
qsub -q highmem Target_Run_AedesL3_merge_noMULE.1k.sh #2 cpu and 40g mem

cd Target
mkdir Target_Run_AedesL3_MULE_msa
cp Target_Run_AedesL3_MULE_finished/AedesL3_MULE_split*/*.msa Target_Run_AedesL3_MULE_msa/
python make_activeTE_sh.py Target/Target_Run_AedesL3_MULE_msa/ '*.msa' activeTE_Target_Run_AedesL3_MULE
python multi_sub_N_move.py . '*_under.sh' 200 aTE_submitted/
python multi_sub_N_move.py . '*_split.sh' 200 aTE_submitted/
python list_aTE_good_bad_mixed.py Target/Target_Run_AedesL3_MULE_msa/ Target/Target_Run_AedesL3_MULE_msa/ AedesL3_MULE_sum

#round2
cat Target/Target_Run_AedesL3_MULE_msa/AedesL3_MULE_sum.conserved_flanks Target/Target_Run_AedesL3_MULE_msa/AedesL3_MULE_sum.bad Target/Target_Run_AedesL3_MULE_msa/AedesL3_MULE_sum.no_tsd > Target/Target_Run_AedesL3_MULE.redo.list
python split_fasta_id.py initial_data_MULE/Merge/AedesL3_MULE.fa
python make_redo_fasta.py Target/Target_Run_AedesL3_MULE.redo.list initial_data_MULE/Merge/AedesL3_MULE.fa
qsub -q highmem Target_Run_AedesL3_MULE_redo.sh #2 cpu and 40g mem

cd Target
mkdir Target_Run_AedesL3_MULE_round2_msa
cp Target_Run_AedesL3_MULE_round2_finished/AedesL3_MULE_split*/*.msa Target_Run_AedesL3_MULE_round2_msa/
python make_activeTE_sh.py Target/Target_Run_AedesL3_MULE_round2_msa/ '*.msa' activeTE_Target_Run_AedesL3_MULE_round2
python multi_sub_N_move.py . '*_under.sh' 200 aTE_submitted/
python multi_sub_N_move.py . '*_split.sh' 200 aTE_submitted/
python list_aTE_good_bad_mixed.py Target/Target_Run_AedesL3_MULE_round2_msa/ Target/Target_Run_AedesL3_MULE_round2_msa/ AedesL3_MULE_round2_sum

