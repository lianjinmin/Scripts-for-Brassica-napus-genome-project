#!/usr/bin/perl
=head1 Name

	auto_ncRNA.pl	-- the pipeline of seach ncRNA.

=head1 Description

	Base on find-RNA.pl and scan_Rfam.pl.
	This program will make directory and work.sh.

=head1 Version

    Version: 1.0    Date: 2010-07-15
    Update: 2016a   Date: 2016-1-20

=head1 Usage

	perl auto_ncRNA.pl [options] <genome.fa>

	note: the <genome.fa> must be input by absolute path

        #rRNA and tRNA
        --rRNA            	search rRNA, default E-value 1e-5
	--evalue              set the E-value for blast in rRNA detection, default E-value 1e-5
        --rnode <str>         set the compute node of rRNA
        --ref_rRNA <file> 	set the reference rRNA seqeunce file, default plant
        --species_type <str> set the type of species to choose the reference rRNA seqeunce file. plant/vertebrate/invertebrate/yeast, default: plant
        --tRNA            	search tRNA
        --tnode <str>         set the compute node of tRNA
        --spec_tag <str>	set non-eukayrote species option for tRNAscan-SE: B,A,O,
        --lines_for_rtRNA <int> set number of lines to form a job, default 1
        
        #miRNA and snRNA
        --miRNA			search miRNA
        --mnode <str>         set the compute node of miRNA
        --snRNA			search snRNA
        --snode <str>         set the compute node of snRNA
        --rfam_dir <str>  	set the directory of Rfam database, default Rfam.12.0
        --lines_for_msRNA <int> set number of lines to form a job, default 100
        --queue <str>     	set the queue in qsub, default no
        --pro_code <str>	set the code for project,default no
        --resource <str>  	set the required resource used in qsub -l option, default vf=1G
        --jobprefix <str> 	set the prefix tag for qsubed jobs, default work
        --step <int>      	set the step number, from where to start, default=1
        --outdir <str>          set the result directory, default="./"
        --prefix <str>    	set a prefix name for the gene ID in gff3 result file
        --run <str>       	set the parallel type, qsub, or multi, default=qsub
        --cpu <int>       	set the cpu number to use in parallel, default=3
        --verbose         	output running progress information to screen
        --help            	output help information to screen

=head1 Exmple

=cut

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use lib "$Bin/../../common_bin";
use GACP qw(parse_config);

my ($rRNA,$ref_rRNA,$species_type,$tRNA,$Spec_tag,$tRNA_seq_file,$lines_for_rtRNA,$evalue);
my ($miRNA,$snRNA,$Rfam_dir,$lines_for_msRNA,$Queue,$Resource,$Jobprefix,$Step,$queue,$Pro_code);
my ($Outdir,$Prefix,$Run,$Cpu,$Verbose,$Help,$Tnode,$Rnode,$Mnode,$Snode);

GetOptions(
	"rRNA"=>\$rRNA,
	"evalue:s"=>\$evalue,
    "species_type:s"=>\$species_type,
	"ref_rRNA:s"=>\$ref_rRNA,
        "rnode:s"=>\$Rnode,
	"tRNA"=>\$tRNA,
#        "trnafile:s"=>\$tRNA_seq_file,
        "tnode:s"=>\$Tnode,
	"spec_tag:s"=>\$Spec_tag,
	"lines_for_rtRNA:i"=>\$lines_for_rtRNA,

	"miRNA"=>\$miRNA,
        "mnode:s"=>\$Mnode,
	"snRNA"=>\$snRNA,
        "snode:s"=>\$Snode,
	"rfam_dir:s"=>\$Rfam_dir,
	"lines_for_msRNA:i"=>\$lines_for_msRNA,
	"queue:s"=>\$queue,
	"pro_code:s"=>\$Pro_code,
	"resource:s"=>\$Resource,
	"jobprefix:s"=>\$Jobprefix,
	"step:i"=>\$Step,

	"outdir:s"=>\$Outdir,
	"prefix:s"=>\$Prefix,
	"run:s"=>\$Run,
	"cpu:s"=>\$Cpu,
	"verbose"=>\$Verbose,
	"help"=>\$Help,
);

die `pod2text $0` if(@ARGV==0 || $Help);

my $config_file = "$Bin/../../config.txt";
my $find_RNA = parse_config($config_file,"find-RNA_path");
my $scan_Rfam = parse_config($config_file,"scan_Rfam_path");
my $stat_prgram = parse_config($config_file,"stat_program_ncRNA");
my $change_to_60 = parse_config($config_file,"change_to_60");

#$ref_rRNA ||= "/ifs2/BC_GAG/Bin/Annotation/bin/Annotation_pipeline1_1.0/04.ncRNA-finding/rRNA-tRNA-Rfam/dat/ref-rRNA/plant-rRNA/plant_rRNA.fa";

## changed by hlj
$species_type ||= "plant";

if($species_type=~/^plant$/i){
   $ref_rRNA ||= parse_config($config_file,"plant_ref_rRNA");
}elsif($species_type=~/^yeast$/i){
   $ref_rRNA ||= parse_config($config_file,"yeast_ref_rRNA");
}elsif($species_type=~/^vertebrate$/i){
   $ref_rRNA ||= parse_config($config_file,"vertebrate_ref_rRNA");
}elsif($species_type=~/^invertebrate$/i){
    $ref_rRNA ||= parse_config($config_file,"invertebrate_ref_rRNA");
}


$lines_for_rtRNA ||= 1;
$evalue ||= 1e-5;
$Rfam_dir ||= parse_config($config_file,"Rfam");
$lines_for_msRNA ||=1000;
#$queue ||= "bc.q";
#$Pro_code ||="gagtest";
$Resource ||="vf=1G";
$Jobprefix ||="work";
my $Rnode_para=(defined $Rnode)?"-node $Rnode":"";
my $Tnode_para=(defined $Tnode)?"-node $Tnode":"";
my $Mnode_para=(defined $Mnode)?"-node $Mnode":"";
my $Snode_para=(defined $Snode)?"-node $Snode":"";
$Step ||= 1;

#$Outdir ||= ".";
$Outdir ||= `pwd`;
chomp($Outdir);
$Run ||= "qsub";
$Cpu ||= 3;

my $seq_file = shift;
my $seq_name = basename($seq_file);
my $species_name=$1 if ($seq_name=~/^(\w+)\./);
#my $tRNA_seq_name = basename($tRNA_seq_file);
my $dir = `pwd`; chomp ($dir);

open (OUT , ">" . "STEP01_ncRNA_finding.sh") or die $!;
print OUT "mkdir ncRNA\ncd ncRNA\n";
if (defined $rRNA)
{
	print OUT "mkdir rRNA\n";
	print OUT "cd rRNA\n";
	#print OUT "nohup perl $find_RNA -rRNA -cpu $Cpu -run $Run -lines $lines_for_rtRNA -outdir $Outdir -ref_rRNA $ref_rRNA -queue $queue ";
	#print OUT "nohup perl $find_RNA -rRNA -cpu $Cpu -run $Run -lines $lines_for_rtRNA -ref_rRNA $ref_rRNA -queue $queue ";
	print OUT "nohup perl $find_RNA -rRNA -cpu $Cpu -run $Run $Rnode_para -lines $lines_for_rtRNA -ref_rRNA $ref_rRNA -evalue $evalue ";
	print OUT "-queue $queue " if (defined $queue);
	print OUT "-pro_code $Pro_code " if (defined $Pro_code);
	print OUT "-spec_tag $Spec_tag " if (defined $Spec_tag);
	print OUT "-prefix $Prefix " if (defined $Prefix);
	print OUT "$seq_file &\n";
	print OUT "cd ..\n";
}
if (defined $tRNA)
{
	print OUT "mkdir tRNA\n";
	print OUT "cd tRNA\n";
	print OUT "perl $change_to_60 $seq_file >$seq_name.60\n";
	#print OUT "nohup perl $find_RNA -tRNA -cpu $Cpu -run $Run -lines $lines_for_rtRNA -outdir $Outdir -queue $queue ";
	#print OUT "nohup perl $find_RNA -tRNA -cpu $Cpu -run $Run -lines $lines_for_rtRNA -queue $queue ";
	print OUT "nohup perl $find_RNA -tRNA -cpu $Cpu -run $Run $Tnode_para -lines $lines_for_rtRNA -evalue $evalue ";
	print OUT "-queue $queue " if (defined $queue);
	print OUT "-pro_code $Pro_code " if (defined $Pro_code);
	print OUT "-spec_tag $Spec_tag " if (defined $Spec_tag);
	print OUT "-prefix $Prefix " if (defined $Prefix);
	print OUT "$seq_name.60 &\n";
	print OUT "cd ..\n";
}
if (defined $miRNA)
{
	print OUT "mkdir miRNA\n";
	print OUT "cd miRNA\n";
	#print OUT "nohup perl $scan_Rfam -type miRNA -cpu $Cpu -run $Run -lines $lines_for_msRNA -outdir $Outdir -queue $queue ";
	#print OUT "nohup perl $scan_Rfam -type miRNA -cpu $Cpu -run $Run -lines $lines_for_msRNA -queue $queue ";
	print OUT "nohup perl $scan_Rfam -type miRNA -cpu $Cpu -run $Run $Mnode_para -lines $lines_for_msRNA -evalue $evalue ";
	print OUT "-resource $Resource -step $Step -jobprefix $Jobprefix -rfam_dir $Rfam_dir ";
	print OUT "-queue $queue " if (defined $queue);
	print OUT "-pro_code $Pro_code " if (defined $Pro_code);
	print OUT "-prefix $Prefix " if (defined $Prefix);
	print OUT "$seq_file &\n";
	print OUT "cd ..\n";
}
if (defined $snRNA)
{
	print OUT "mkdir snRNA\n";
	print OUT "cd snRNA\n";
	#print OUT "nohup perl $scan_Rfam -type snRNA -cpu $Cpu -run $Run -lines $lines_for_msRNA -outdir $Outdir -queue $queue ";
	#print OUT "nohup perl $scan_Rfam -type snRNA -cpu $Cpu -run $Run -lines $lines_for_msRNA -queue $queue ";
	print OUT "nohup perl $scan_Rfam -type snRNA -cpu $Cpu -run $Run $Snode_para -lines $lines_for_msRNA -evalue $evalue ";
	print OUT "-resource $Resource -step $Step -jobprefix $Jobprefix -rfam_dir $Rfam_dir ";
	print OUT "-queue $queue " if (defined $queue);
	print OUT "-pro_code $Pro_code " if (defined $Pro_code);
	print OUT "-prefix $Prefix " if (defined $Prefix);
	print OUT "$seq_file &\n";
	print OUT "cd ..\n";
}
print OUT "cd ..\n";
close (OUT);

open (OUT , ">" . "STEP02_ncRNA_statistics.sh") or die $!;
print OUT "mkdir ncRNA_statistics\ncd ncRNA_statistics\n";
my $stat_parameter = "";
if (defined $miRNA)
{
	print OUT "ln -s $dir/ncRNA/miRNA/$species_name.miRNA.gff miRNA.gff\n";
	$stat_parameter = $stat_parameter . "-miRNA" . " ";
}
if (defined $tRNA)
{
	print OUT "ln -s $dir/ncRNA/tRNA/$species_name.tRNA.gff tRNA.gff\n";
	$stat_parameter = $stat_parameter . "-tRNA" . " ";
}
if (defined $rRNA)
{
	print OUT "ln -s $dir/ncRNA/rRNA/$species_name.rRNA.gff rRNA.gff\n";
	$stat_parameter = $stat_parameter . "-rRNA" . " ";
}
if (defined $snRNA)
{
	print OUT "ln -s $dir/ncRNA/snRNA/$species_name.snRNA.gff snRNA.gff\n";
	$stat_parameter = $stat_parameter . "-snRNA" . " ";
}
print OUT "ln -s $seq_file $seq_name\n";
print OUT "perl $stat_prgram $stat_parameter $seq_name\n";
close (OUT);

open (OUT , ">" . "STEP03_delete_tmp_file.sh") or die $!;
	print OUT "cd ncRNA\n";
	if (defined $rRNA)
	{
		print OUT "cd rRNA\n";
		print OUT "rm -r *.cut  *.shell *.shell.*\n";
		print OUT "cd ..\n";	
	}
	if (defined $tRNA)
	{
		print OUT "cd tRNA\n";
		print OUT "rm -r *.cut *.shell *.shell.*\n";
		print OUT "cd ..\n";
	}
	if (defined $miRNA)
	{
		print OUT "cd miRNA\n";
		print OUT "rm -r *.shell *.shell.* Rfam.fasta.* *.cmsearch\n";
		print OUT "cd ..\n";	
	}
	if (defined $snRNA)
	{
		print OUT "cd snRNA\n";
		print OUT "rm -r *.shell *.shell.* Rfam.fasta.* *.cmsearch\n";
		print OUT "cd ..\n";
	}
	print OUT "cd ..";
close (OUT);
