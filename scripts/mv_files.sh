mkdir bkp
for i in {26..50}; do
  arquivo="S$i"  
  mv $(ls| grep $arquivo) bkp/ 
done