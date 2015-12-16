#*** Initialization ***
my $primer='';
my $Filename='';
my @file_data=();
my $pmid='';
my $filename='';
my @file=();

#*** Initialization ***
   my @sentences=();
    my $line='';
    my $temp='';
    my $empty='';
    my $sentence='';
    my $pmidSeq='';
    
#*** Initialization ***
   my @readblast=();
   my $line2='';
    my $blastfile='';

#*** Open input file***
$filename='out1.txt';
open(FILE,$filename)or die("Couldn't open Input file\n");
@file_data = <FILE>;
close(FILE);
$pmidSeq='>>out2.txt';
open(OUTPUT, $pmidSeq)or die("Couldn't open Output file\n");;	

 # print  @file_data;
foreach $line (@file_data)
     {
        if ($line =~ /^     \d/){
          
          @sentences = split(':', $line);
        # *** searching for primer in each small fragment ***
	     $temp = @sentences;
              print OUTPUT $temp ;
          }
}
           
            close(OUTPUT);
  exit;
