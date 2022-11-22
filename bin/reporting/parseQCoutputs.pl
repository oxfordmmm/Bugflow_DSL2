#!/usr/bin/perl

#use PDF::API2;
#use DateTime;
use Time::Piece;

$fqc = $ARGV[0]; #FastQC output
$qst = $ARGV[1]; #Quast output

open(qc, $fqc);
open(qt, $qst);

$outqc = "QC_summary_table.tsv";
$sample = "sample_name.tsv";

system("echo Sample Name\t$fqc | cut -f 1 -d '.' ");
#system("echo Sample Name\t$fqc | cut -f 1 -d '.' >> $outqc");

open(outqc, '>>', $outqc);
open(sample, '>', $sample);

#print outqc "Sequence Quality\n";
print "\nGENOME QUALITY STATS\n";
print "\nMetric\tExpected value\tActual value\tQuality\n\n";
print outqc "\nMetric\t\t\tExpected value\t\t\tActual value\t\t\tQuality\n";

while(<qc>){
    chomp;
    $f = $_;
    #print "$f\n";
    #if ($f =~ "Basic Statistics" && $f =~ "pass"){print "Overall sequence quality\tPass\n"} elsif ($f =~ "Basic Statistics" && $f =~ "fail") {print "Overall sequencing quality\tFail\n"}
    if ($f =~ /Filename/){
        $date = localtime->strftime('%m/%d/%Y'); 
        #print "$datetime\n";
        @fname = split/\t/,$f;
        for $n (0..$#fname){}
        $fname[1] =~ s/\_raw\_reads\.fastq\.gz//g;
        #print sample "Sample details\n";
        print sample "Specimen identifier:\t\t$fname[1]\t\t\t\tCollection date:\t\t$date\n";
        print sample "Sequence GUUID:\t\t$fname[1]\t\tSequencing date:\t\t$date\n";
        print sample "Sequence plate name:\t\t$fname[1]\t\t\t\tReport date:\t\t$date\n";
    }
    if ($f =~ /Total Sequences/){
        @seqs = split/\t/,$f;
        for $s (0..$#seqs){}
        if ($seqs[1] == 1000000 && $seqs[1] < 20000000 && $seqs[1] == 20000000){
            print "Total Sequences\t1.0-20.0 M\t$seqs[1]\tPass\n";
            print outqc "Total Sequences\t\t\t1.0-20.0 M\t\t\t$seqs[1]\t\t\tPass\n";
        }
        else {
            print "Total Sequences\t1.0-20.0 M\t$seqs[1]\tFail\n";
            print outqc "Total Sequences\t\t\t1.0-20.0 M\t\t\t$seqs[1]\t\t\tFail\n";
        }
    }
    if ($f =~ /\%GC/){
        @gc = split/\t/,$f;
        for $g (0..$#gc){}
        if ($gc[1] == 27.9 && $gc[1] < 29.2 && $gc[1] == 29.2){
            print "%GC\t27.9-29.20\t$gc[1]\tPass\n";
            print outqc "%GC\t\t\t27.9-29.20\t\t\t$gc[1]\t\t\tPass\n";
        }
        else {
            print "%GC\t27.9-29.20\t$gc[1]\tFail\n";
            print outqc "%GC\t\t\t27.9-29.20\t\t\t$gc[1]\t\t\tFail\n";
        }
    }

    if ($f =~ /Sequence length/){
        @len = split/\t/,$f;
        for $s (0..$#len){}
        #print "$len[1]*\n";
        @seqlen = split/-/,$len[1];
        for $s (0..$#seqlen){}
        $mean = ($seqlen[0] + $seqlen[1])/2;
        #print "$mean**\n";
        if ($mean > 50 | $len[1] > 150){
            print "Sequence length (bp)\t50-150M\t$mean\tPass\n";
            print outqc "Sequence length (bp)\t\t\t50-150M\t\t\t$mean\t\t\tPass\n";
        }
        else {
            print "Sequence length (bp)\t50-150M\t$len[1]\tFail\n";
            print outqc "Sequence length (bp)\t\t\t50-150M\t\t\t$len[1]\t\t\tFail\n";
        }
    }

}

close(qc);

while(<qt>){
    chomp;
    $q = $_;
    #print "$q\n";
    if ($q =~ "Total length" && $q =~ "(>= 0 bp)"){
        @gnom = split/\t/,$q;
        for $m (0..$#gnom){}
        if ($gnom[1] == 3900000 && $gnom[1] < 4500000 && $gnom[1] == 4500000){
            print "Total assembly size (bp)\t3.9-4.5\t$gnom[1]\tPass\n";
            print outqc "Total assembly size (bp)\t\t\t3.9-4.5\t\t\t$gnom[1]\t\t\tPass\n";
        }
        else {
            print "Total assembly size (bp)\t3.9-4.5\t$gnom[1]\tFail\n";
            print outqc "Total assembly size (bp)\t\t\t3.9-4.5\t\t\t$gnom[1]\t\t\tFail\n";
        }
    }
    if ($q =~ "Largest contig"){
        @con = split/\t/,$q;
        for $c (0..$#con){}
        print  "Largest contig (bp)\tNA\t$con[1]\tNA\n";
        print outqc  "Largest contig (bp)\t\t\tNA\t\t\t$con[1]\t\t\tNA\n";
    }
    if ($q =~ "N50"){
        @nfifty = split/\t/,$q;
        for $s (0..$#nfifty){}
        print "N50 (bp)\tNA\t$nfifty[1]\tNA\n";
        print outqc "N50 (bp)\t\t\tNA\t\t\t$nfifty[1]\t\t\tNA\n";
    }
    if ($q =~ /\%GC/){
        @gc = split/\t/,$q;
        for $g (0..$#gc){}
        if ($gc[1] == 27.9 && $gc[1] < 29.2 && $gc[1] == 29.2){
            print outqc "%GC\t\t\t27.9-29.20\t\t\t$gc[1]\t\t\tPass\n";
        }
        else {
            print outqc "%GC\t\t\t27.9-29.20\t\t\t$gc[1]\t\t\tFail\n";
        }
    }


}
print "\n\n";

close(sample);

close(qt);

close(outqc);


exit;
