for t in luad breast skin colon  kidney prostate; do #ovary
   echo ${t}
   if [ ! -d data/procData/${t}/methylation/bins/ ]
      then
         mkdir data/procData/${t}/methylation/bins/
   fi
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
   for bin in bins10kb bins100kb bins1Mb ; do
      echo ${bin}
      bedmap --echo --wmean --delim '\t' \
         data/procData/${bin}_allChrs.bed \
         data/procData/${t}/methylation_methbank/${T}_liftover_sorted.bed \
         > data/procData/${t}/methylation/bins/${bin}_methylation.out.bed
   done  
done
