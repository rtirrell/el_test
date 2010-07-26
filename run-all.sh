for f in `find data -name genome_*`; do 
  echo $f && python el_test.py $f ceu; 
done

