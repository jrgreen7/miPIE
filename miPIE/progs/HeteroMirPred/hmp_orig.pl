#!/usr/bin/perl -w

use strict;
use warnings;

#*** Initialization ***
   my @file_data=();
   my $filename='';
   my @file=();
   my @sentences=();
    my $line='';
    my $tempfile='';
    my $temp='';
    my $tempdata='';
    my $tempdata2='';
    my $empty='';   
    my $sentence='';
    my $fastaSeq='';
   my @readspec=();
   my $foldfile='';
   my $allstatfile='';
   my $num = 0 ;

my $statfilename='';
my $feature20name='';
my $SCfilename='';
my $spectralfile='';
my $mfe14file = '';
my $mfe23file ='';
my $BP1file ='';
my $BP2file ='';
my $line1='';
my $line2='';
my $line3='';
my $line4='';
my $line5='';
my $line6='';
my $line7='';
my $line8='';
my $line9='';
my $spec='';
my $tab = "\t";


#*** Open input file***

$filename =  $ARGV[0];

open(FILE,$filename)or die("Couldn't open Input file\n");
@file_data = <FILE>;
close(FILE);
system("rm $filename.csv");
#system("rm $filename.20features.csv");

open(FILE3, ">> $filename.csv");
$tempdata = "mfe,Prob,efe,mfe1,mfe2,dG,dP,dQ,dF,zG,zP,zQ,mfe3,mfe4,nefe,freq,div,diff,dH,dHL,dH/Loop,dS,dSL,dS/loop,Tm,TmL,Tm/Loop,au/L,gc/L,gu/L,P_au/bp,P_gc/bp,P_gu/bp,tot_bp/Loops,mfe5,SCI,SC/non_tot,SCxMFE/Mean,SCxabsZG,SCxdP,SC/len,SC/auLgcLguL,SC/nonauLgcLguL,SC/1dP,SC/nonA,SC/nonU,score1,score2,score3,score4,score5,score6,score7,score8,score9,score10,score11,score12,score13,score14,score15,score16,score17,score18,score19,score20,score21,score22,score23,score24,score25,score26,score27,score28,score29,score30,score31,score32,score33,score34,score35,score36,score37,score38,score39,score40,score41,score42,score43,score44,score45,score46,score47,score48,score49,score50,score51,score52,score53,score54,score55,score56,score57,score58,score59,score60,score61,score62,score63,score64,score65,score66,score67,score68,score69,score70,score71,score72,score73,score74,score75,score76,score77,score78,score79,score80,score81,score82,score83,score84,score85,score86,score87,score88,score89,score90,score91,score92,score93,score94,score95,score96,score97,score98,score99,score100,a1,a2,a3,a4,a5,a6,a7,a8,g1,g2,g3,g4,g5,g6,g7,g8,c1,c2,c3,c4,c5,c6,c7,c8,t1,t2,t3,t4,t5,t6,t7,t8,mfe_a,MeanBP,A,C,G,U,AA,AC,AG,AU,CA,CC,CG,CU,GA,GC,GG,GU,UA,UC,UG,UU,pairprob1,pairprob2,pairprob3,pairprob4,pairprob5,pairprob6,pairprob7,pairprob8,pairprob9,pairprob10,Non_BPP,nonBP_a,nonBP_c,nonBP_g,nonBP_u,Class\n";
print FILE3 $tempdata;

#open(FILE4, ">> $filename.20features.csv");
#$tempdata2 = "Prob,mfe1,zG,zP,dH/Loop,Tm/Loop,au/L,P_gu/bp,tot_bp/Loops,mfe5,SCI,scxdP,SCxMFE/Mean,SC/1dP,SC/non_A,Non_BPP,a8,c8,g3,pairprob4\n";
#print FILE4 $tempdata2;

 # print  @file_data;
foreach $line (@file_data)
     {
        if ($line =~ /^>/){
          $temp  .= $line;
          $num++;
	  $fastaSeq=">$filename.fasta";
	  #$foldfile=">>$num.fold";	                   
          }
       
       elsif ($line =~ /^\w/){
             $temp .= $line;
          
	     open(OUTPUT, $fastaSeq);	
             print OUTPUT $temp;
	     
   	     print("\n-calculating features-------------\n");

	     system("RNAfold < $filename.fasta > $filename.fold"); 
             system("perl genRNAStats_mod.pl < $filename.fasta > $filename.stat2");
	     system("perl genRNAStats.pl < $filename.fasta > $filename.data1");
	     system("./RNAtopological.exe < $filename.fold > $filename.spec");
	     
	     system("perl genRandomRNA.pl -n 50 -m d < $filename.fasta > $filename.random.fasta");
		system("RNAfold < $filename.random.fasta > $filename.random.fold"); 
		system("rm *_ss.ps");
	     
 	     system("perl genRNARandomStats_mod.pl -n 50 -i $filename.random.fold -o $filename.zdata -m $filename.fold");
	     $tempfile=">temp.txt";
		system("rm *.random.fold");
   		system("rm *.random.fasta");
   		system("java mfe14 $filename.data1");
   		system("java mfe23 $filename.spec");
		system("perl melt2.pl $filename.fasta > $filename.mfold2");
		system("perl melt.pl $filename.fasta > $filename.mfold");
		system("java Mfoldfilter $filename");
		system("python selfcontain.py -i $filename.fasta > $filename.sc");
   		system("java bpcount1_mod $filename");
                #Count Percent of Trimer occurs in LongStem------------------------------- 
                #system("perl trimerstemcount_2.pl < $filename.longstem > $filename.trimer");
   		system("RNAfold -p2 < $filename.fasta > $filename.RNAfold1");  
   		system("java RNAfoldfilter $filename");
   		system("rm *.ps");
		system("rm $filename.fasta.37.ext");
		system("rm $filename.fasta.ct");
		system("rm $filename.fasta.run");
		system("./randfold -d $filename.fasta 100 > $filename.prob");
                system ("./Conv.exe $filename.fasta $filename.f2 1");
		system ("./Extract_PP.exe -f $filename.f2 -out $filename.Probpair");
		system("java filter $filename"); 
		print(" Features for sequence no. $num is Done ----------------------------\n");	
		
		my $statfilename= "$filename.feature";
		
		open(FILE1, $statfilename);
		
                #my $feature20name= "20_features.txt";

		while(defined(my $line1 = <FILE1>)) {
			print FILE3 $line1;
			}
                #open(FILE2, $feature20name);
		#$line2 = <FILE2>;	
		#print FILE4 $line2;

	         system("rm $filename.fasta");
                 system("rm $filename.feature");
		 system("rm $filename.bp1");
 		 system("rm $filename.data1");
		 system("rm $filename.longstem");
		 system("rm $filename.mfold2");
 		 system("rm $filename.mfold");
                  system("rm $filename.sc");
		 system("rm $filename.RNAfold1");
		 system("rm $filename.RNAfold2");
		 system("rm $filename.spec");
		 system("rm $filename.stat2");
		 system("rm $filename.fasta.37.plot");
		 system("rm $filename.spec.mfe23");
		 system("rm $filename.data1.mfe14");	
		system("rm $filename.fasta.dG");
                 system("rm $filename.fold");	
                system("rm $filename.zdata");
		system("rm $filename.prob");
                system("rm $filename.f2");
                system("rm $filename.Probpair");
		#my $feature20name="$num.predict";
		
           #  print TEMPOUT $num;
             close(OUTPUT);
            # close(TEMPOUT);
             $temp=$empty;

           }
	

}
$tempdata="-19.90,0.712871,-21.52,-0.0043260422508880165,-0.04627906976744185,-0.2314,0.2674,0.4940,0.267949,0.5555,-0.0085,0.1494,-0.038565891472868215,-0.8652173913043477,-0.2502325581395349,0.0725996,13.76,0.018837209302325592,-219.6,-2.5534883720930233,-36.6,-649.7,-7.554651162790698,-108.28333333333335,64.9,0.7546511627906978,10.816666666666668,0.08139534883720931,0.13953488372093023,0.046511627906976744,30.434782608695652,52.17391304347826,17.391304347826086,3.8333333333333335,-0.01617531111056922,0.233076923077,0.0036996336996349207,0.2144767254498007,0.1294742307692735,0.062324769230789806,0.0027101967799651166,0.8715050167226955,0.3181684981686032,0.5010252000795358,0.4905767174558206,0.40488221134651,31.689809,32.621971,32.362093,41.740111,31.749816,33.145429,31.654241,33.319215,34.239531,33.798598,33.063755,32.584120,30.743285,29.428007,30.112283,37.443997,39.103310,32.459125,34.525427,40.814362,29.224489,29.396110,31.437623,39.257580,25.080382,31.601672,35.393598,39.371427,33.840741,33.948327,36.371919,34.266313,34.810444,34.591607,34.576406,39.009769,34.489710,33.972398,36.510383,32.770041,35.886874,32.132604,29.695694,36.074119,32.736164,26.738074,35.754237,37.456652,33.470867,32.835982,28.862284,37.810337,32.248926,34.439803,44.317952,37.187071,32.218946,39.013395,43.445849,31.949933,38.853779,32.553398,32.643886,39.741112,36.927307,43.484368,36.506700,37.681005,42.601370,39.826091,36.686535,35.365687,39.809463,36.274141,37.024240,42.756308,35.397101,34.322161,35.354594,37.010754,42.539842,38.063533,35.511423,45.156984,40.640390,47.962703,34.986770,38.522193,37.625383,32.706787,33.048025,35.316596,38.739237,36.507305,35.132372,35.916107,35.664156,36.276664,36.516215,37.187273,0.144320,0.077597,0.092491,0.160699,0.060964,0.054579,0.117236,0.292114,0.198193,0.175294,0.059605,0.095407,0.080521,0.092747,0.084924,0.213308,0.048891,0.060907,0.038149,0.054893,0.197528,0.099743,0.207573,0.292317,0.185159,0.106363,0.061837,0.196743,0.075486,0.065256,0.097120,0.212035,-24.950000,26.222706,0.162791,0.302326,0.232558,0.302326,0.023529,0.047059,0.035294,0.058824,0.082353,0.094118,0.023529,0.105882,0.035294,0.058824,0.058824,0.082353,0.023529,0.105882,0.117647,0.047059,1.898788,2.908062,1.852185,2.616498,2.330157,2.652784,3.893344,2.026430,5.318997,1.827623,0.458332,0.475108,0.528499,0.202839,0.575666,Other";
print FILE3 $tempdata;

system("java -classpath \$CLASSPATH:weka.jar:libsvm.jar:mysql-connector-java-5.0.8-bin.jar weka.classifiers.meta.Vote -l Vote.model -T $filename.csv -p 36 > out1.txt");

system("clear");

print("\n\nInput file: $filename; There are $num query sequences in the input file. \n");
print("\n.... Completed Extracting All Features ....\n");
print("\nAll features extracted are written to the file: $filename.csv.");

print("\n\n############ Prediction results by Heterogeneous Ensemble #############\n\n");


system("grep 1 out1.txt | awk '{gsub(/+/,\" \")}; {print \"  Sequence no. \"\$1 \" predict as \" \$3  \" - with a probability of \"\$4 \n}' |head -n -1 ");

exit;
