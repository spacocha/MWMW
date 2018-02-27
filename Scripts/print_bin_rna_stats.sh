# print tRNAs
for i in `ls *.gff`; do 
    echo $i; 
    cat $i | grep "Aragorn:1.2" | python -c 'from sys import stdin, stdout; d = stdin.read().split("\n")[:-1]; dd = [i.split("product=")[-1] for i in d]; print len(set(dd))'; 
done

for i in `ls *.gff`; do
    echo $i;
    cat $i | grep "23S ribosomal RNA" | wc -l;
    cat $i | grep "16S ribosomal RNA" | wc -l;
    cat $i | grep "5S ribosomal RNA" | wc -l;
done

