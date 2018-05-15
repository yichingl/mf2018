#/bin/csh -f

set tmpf = /tmp/acc.$$
set tmpf2 = /tmp/acc2.$$
onintr cleanup

cat /dev/null >! $tmpf2
foreach x (`interval 1 10`)
    set setout = SET.$x.RESULTS
    awk '{print "1 ",$0}' SET.$x.class1 >! $tmpf
    awk '{print "0 ",$0}' SET.$x.class0 >> $tmpf
    cat $tmpf >> $tmpf2
    awk '{if($1==1){N1++;N++;if($4==1){N1c++;Nc++;}}if($1==0){N0++;N++;if($4==0){N0c++;Nc++}}}END{printf("recogn->\tclass0\tclass1\nclass0\t%d\t%d\nclass1\t%d\t%d\n\nClass 0 accur = %f (%d/%d)\nClass 1 accur = %f (%d/%d)\nOverall accur = %f (%d/%d)\n",N0c,(N0-N0c),(N1-N1c),N1c,100*N0c/N0,N0c,N0,100*N1c/N1,N1c,N1,100*Nc/N,Nc,N)}' $tmpf >! $setout
end

awk '{if($1==1){N1++;N++;if($4==1){N1c++;Nc++;}}if($1==0){N0++;N++;if($4==0){N0c++;Nc++}}}END{printf("recogn->\tclass0\tclass1\nclass0\t%d\t%d\nclass1\t%d\t%d\n\nClass 0 accur = %f (%d/%d)\nClass 1 accur = %f (%d/%d)\nOverall accur = %f (%d/%d)\n",N0c,(N0-N0c),(N1-N1c),N1c,100*N0c/N0,N0c,N0,100*N1c/N1,N1c,N1,100*Nc/N,Nc,N)}' $tmpf2 >! OVERALL.RESULTS

cleanup:
/bin/rm -f $tmpf $tmpf2

exit -0
