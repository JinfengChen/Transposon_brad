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

echo "Run target"
#python make_target_nonauto_general.py /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/NonAutonomous/initial_data/AaegL3.MH_RM/query /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/NonAutonomous/Target /rhome/cjinfeng/BigData/00.RD/Mosquito_TE/reference/AaegL3.fa Target_Run
#qsub Target_Run_AedesL3_unique.sh
#qsub Target_Run_AedesL3_bad-to-unique_unique.sh

echo "Run activeTE"
mkdir Target/Group_msa
cp Target/Target_Run_2015_04_30_130853/AedesL3_unique_split*/*.msa Target/Group_msa
python make_activeTE_sh.py Target/Group_msa/ '*.msa' activeTE_test

