rm tmp
touch tmp
for i in `seq 0 50 3100`; do cat tmp D2O32_md_polar-HOMO_centers_s1-1_$i.xyz > tmp0; mv tmp0 tmp; done
mv tmp HOMO_centers.xyz
