for f in `find data/ -name "genome_*"`; do 
  echo $f && python el.py $f c; 
done

