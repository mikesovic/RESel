cd ${1}_out
mkdir feature_files
mkdir summary_files

touch features_detail_ALL_head.txt
echo -e 'feature_ID\tfeature_length\tsites_covered\tpct_covered\tRE1\tRE2\toverhang1\toverhang2' >> features_detail_ALL_head.txt

cat   *features_detail.txt > features_tmp1.txt
mv *features_detail.txt feature_files

grep -v -E "overhang2$" features_tmp1.txt > features_tmp2.txt
cat features_detail_ALL_head.txt features_tmp2.txt > all_features_detail.txt


touch summary_ALL_head.txt
echo -e 'genome\tRE1\tRE2\tadaptor1\tadaptor2\tnum_in_range\tfeatures_evaluated\tfeatures_covered_partial\tfeatures_covered_whole\ttotal_features_covered\tpct_features_covered' >> summary_ALL_head.txt

cat *summary.txt > summary_tmp1.txt
mv *summary.txt summary_files

grep -v -E "pct_features_covered$" summary_tmp1.txt > summary_tmp2.txt
cat summary_ALL_head.txt summary_tmp2.txt > all_summary.txt

rm features_tmp*.txt
rm summary_tmp*.txt

#for i in `ls *features_detail.txt`
#do
#tail -n +2 ${i} >> features_detail_ALL.txt 
#done

#for i in `ls *summary.txt`
#do
#tail -n +2 ${i} >> all_summary.txt
#done

#mv *features_detail.txt feature_files
#mv *summary.txt summary_files
#mv feature_files/all_features_detail.txt ./
#mv summary_files/all_summary.txt ./

