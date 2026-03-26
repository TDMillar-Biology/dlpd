while read line
  do strain=$(echo $line | cut -f1 -d' ')
  echo -e "contig\tbreak" > ${strain}.tsv
done < ../master_control.tsv
