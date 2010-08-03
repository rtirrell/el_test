for f in `find data -name genome_*.txt`; do 
  echo $f && python longevity.py $f ceu codename-`basename $f`; 
done

