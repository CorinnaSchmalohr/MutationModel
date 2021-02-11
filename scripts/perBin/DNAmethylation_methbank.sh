for t in luad breast skin colon  kidney prostate; do #ovary
   echo ${t}
   case $t in
   luad)
   T=Lung
   ;;
   breast)
   T=Breast
   ;;
   skin)
   T=Skin
   ;;
   colon)
   T=Colon
   ;;
   kidney)
   T=Kidney
   ;;
   prostate)
   T=Prostate
   ;;
   esac
   echo $T
   for bin in bins10bp bins100bp bins1kb bins10kb bins100kb bins1Mb ; do
      echo ${bin}
      mkdir data/procData/bins/${bin}/${t}/methylation/
      bedmap --echo --wmean --delim '\t' \
         data/procData/bins/${bin}/${bin}.bed \
         data/procData/${t}/methylation_methbank/${T}_liftover_sorted.bed \
         > data/procData/bins/${bin}/${t}/methylation.out.bed
      bedmap --echo --wmean --delim '\t' \
         data/procData/bins/${bin}/${bin}_${t}_covered.bed \
         data/procData/${t}/methylation_methbank/${T}_liftover_sorted.bed  \
         > data/procData/bins/${bin}/${t}/methylation_covered.out.bed 
   done  
done
