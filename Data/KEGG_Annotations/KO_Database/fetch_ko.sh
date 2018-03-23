while read TO_num; do
    echo $TO_num;
    wget http://rest.kegg.jp/link/ko/$TO_num;
    mv $TO_num TO_Files/$TO_num.txt;
    sleep 3;
done < ./T0_Numbers.txt
