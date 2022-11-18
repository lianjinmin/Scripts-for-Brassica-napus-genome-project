#!/usr/bin/perl

=head1 Name

=head1 Description

=head1 Usage

    perl auto_gene_finding.pl [options] <genome.fa> <remask_genome.fa>
    Options:
        --Homolog
            --Hcpu <i>          set the cpu number to use in parallel, default=100
            --Hrun <s>          set the parallel type, qsub, or multi, default=qsub
            --Htophit <i>       select best hit for each protein, default no limit
            --blast_eval <s>    set eval for blast alignment, default 1e-5
            --align_rate <f>    set the aligned rate for solar result, default 0.25
            --extend_len <i>    set the extend length for genewise DNA fragment, default 500
            --genblasta_opt <s> set the genblasta options, default Genblasta_opt="-p T -e 1e-2 -g T -f F -a 0.5 -d 100000 -r 100 -c 0.5 -s -100 "
            --filter_rate <f>   set the filter rate for the best hit of geneblastA result,default 0.7
            --net <s>           set the net directory result.
            --rgene <s>         set the reference gene file with gff.
            --step <s>          set which step to run, default 1234
            --lines <i>         set the --lines option of qsub-sge.pl,required.must define
            --species <s>       set the homolog species
            --inputpep <s>      set the homolog species protein sequences
            --Hresource <s>     set the qsub-sge commond.default vf=1.0G
            --genblasta         set if use genblasta instead of solar, defaut no.
            --Hremask           select if use the remasked genome
        --Denovo
            --augustus                  run augustus
            --species_for_augustus <s>  set augustus specises
            --genescan                  run genescan
            --genescan_para <s>         give genscan parameter file
            --glimmerhmm                run glimmerhmm
            --glimmerhmm_para <s>       give glimmerHMM parameter file
            --snap                      run snap
            --snap_para <s>             give snap parameter file
            --fgenesh                   run fgenesh
            --geneid                    run geneid
            --geneid_para <s>           give geneid parameter file
            --prefix_for_Denovo <s>     prefix for Denovo
            --Dcpu <s>                  set the cpu number to use in parallel, default=20
            --Drun <s>                  set the parallel type, qsub, or multi, default=qsub
        --EST
            --EST_species <s>   set the EST species
            --input_EST <s>     set the EST species fasta file
            --Ecpu <i>          set the cpu number to use in parallel, default=100
            --Erun <s>          set the parallel type, qsub, or multi, default=qsub
            --identity <s>      set identity cutoff, default 0.95
            --alignrate <s>     set alignrate cutoff, default 0.95
            --dbcut <i>         set the number to cut database file, default=3
            --cluster <s>       run cluster before assembly,T/F, default: T
            --tophit <i>        set the number of top hits for the result, default no limitation
        --cDNA
            --cDNA_species <s>  set the cDNA species
            --input_cDNA <s>    set the EST species fasta file
        --Glean
            --genome_name:s
            --parameter_yaml:s
            --minlen:i
            --minintron:i
            --maxintron:i
            --cds:s
            --scaf:s
            --resource_glean:s
            --cpu_for_glean:i
        --Process
            --SP_cuts:i
            --SP_cpu:i
            --blast_evalue:s
            --blast_vf:s
            --rna_list:s
            --soap_vf:s
            --Cover_cutoff:i
        --CCG
            --reads_list:s
            --train_gff:s
            --combine_gff:s
	    --qual_system:s
            --alignment:s
            --build_vf:s
            --align_vf:s
            --max_i_len:i
            --asmb_vf:s
            --asmb_proc_num:i
            --min_cds_len:i
            --orf:s
            --CCG_step:s
            --CCG_name:s
        --Verbose
        --Help
        --queue:s
        --pro_code:s

=head1 Exmple

	perl /ifs2/BC_GAG/Bin/Annotation/bin/Annotation_pipeline1_1.0/Auto_pipline/auto_gene_finding/bin/auto_geno_finding.pl -queue bc.q -pro_code gagtest -Homolog -species bee -inputpep /panfs/GAG/liushiping/niLztvD/data/homolo/bee.pep.pure.30.fa -species bee2 -inputpep /panfs/GAG/liushiping/niLztvD/data/homolo/bee.pep.pure.30.fa /panfs/GAG/liushiping/niLztvD/data/niLztvD.scafSeq_fill.fa

	perl /ifs2/BC_GAG/Bin/Annotation/bin/Annotation_pipeline1_1.0/Auto_pipline/auto_gene_finding/bin/auto_geno_finding.pl -queue bc.q -pro_code gagtest -Denovo -augustus -species_for_augustus fly -species_for_augustus NL_lsp -genescan /nas/GAG_02/liushiping/GACP/GACP-7.0/01.gene_finding/denovo-predict/dat/genscan_para/HumanIso.smat -glimmerhmm /nas/GAG_02/liushiping/GACP/software/GlimmerHMM/trained_dir/human -snap /panfs/GAG/liushiping/niLztvD/gene/train/snap/NL.hmm /panfs/GAG/liushiping/niLztvD/data/niLztvD.scafSeq_fill.fa

	perl /ifs2/BC_GAG/Bin/Annotation/bin/Annotation_pipeline1_1.0/Auto_pipline/auto_gene_finding/bin/auto_geno_finding.pl -queue bc.q -pro_code gagtest -EST -prefix_for_EST NL -EST_data /panfs/GAG/liushiping/niLztvD/data/EST/Nilaparvata_EST.fa /panfs/GAG/liushiping/niLztvD/data/niLztvD.scafSeq_fill.fa

=cut

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use lib "$Bin/../../common_bin";
use GACP qw(parse_config);


my ($Homolog,$Denovo,$EST,$cDNA,$Glean,$Process,$CCG);
my ($Hcpu,$Hrun,$Hresource,$Htophit,$Blast_eval,$Align_rate,$Extend_len,$Genblasta_opt,$Filter_rate,$Net,$Ref_gene,$Step,$Lines,@Species,@Inputpep,$Genblasta,$Hremask,$i); ## for Homolog
my ($augustus,$genescan,$glimmerhmm,$snap,$fgenesh,$geneid); ## for Denovo
my ($species_for_augustus,$genescan_para,$glimmerhmm_para,$snap_para,$fgenesh_para,$geneid_para,);
my ($prefix_for_Denovo,$Dcpu,$Drun,$Shape,$Denovo_parameter);
my ($Ecpu,$Erun,$DBcut,$Identity_cutoff,$Alignrate_cutoff,$Prefix_for_EST,$Best_hit,$Cluster,@EST_species,@input_EST,@input_EST_name,$EST_parameter,@cDNA_species,@input_cDNA,@input_cDNA_name,$Prefix_for_cDNA); ## for EST
my ($genome_name,$all_glean_gff,$parameter_yaml,$run_for_glean,$minlen,$minintron,$maxintron,$cpu_for_glean,$cds,$scaf,$Resource_glean,$evidence); ## for glean
my ($SP_cuts, $SP_cpu, $blast_evalue, $blast_vf, $rna_list, $soap_vf, $Cover_cutoff); ## for Process
my ($reads_list,$train_gff,$gff_to_combine,$tscore,$tstart,$tstop,$qual_system,$alignment,$build_vf,$align_vf,$max_i_len,$asmb_vf,$asmb_proc_num,$min_cds_len,$CCG_orf,$CCG_step,$CCG_name);  ## for CCG
my ($Verbose,$Help,$all_stat);
my ($Queue,$Pro_code,$opt);

GetOptions(
    "Homolog"=>\$Homolog,
        "Hcpu:i"=>\$Hcpu,
        "Hrun:s"=>\$Hrun,
        "Htophit:i"=>\$Htophit,
        "blast_eval:s"=>\$Blast_eval,
        "align_rate:f"=>\$Align_rate,
        "extend_len:i"=>\$Extend_len,
        "genblasta_opt:s"=>\$Genblasta_opt,
        "filter_rate:f"=>\$Filter_rate,
        "net:s"=>\$Net,
        "rgene:s"=>\$Ref_gene,
        "step:s"=>\$Step,
        "lines:i"=>\$Lines,
        "species:s"=>\@Species,
        "inputpep:s"=>\@Inputpep,
        "Hresource:s"=>\$Hresource,
        "genblasta"=>\$Genblasta,
        "Hremask"=>\$Hremask,
    "Denovo"=>\$Denovo,
        "augustus"=>\$augustus,
        "species_for_augustus:s"=>\$species_for_augustus,
        "genescan"=>\$genescan,
        "genescan_para:s"=>\$genescan_para,
        "glimmerhmm"=>\$glimmerhmm,
        "glimmerhmm_para:s"=>\$glimmerhmm_para, 
        "snap"=>\$snap,
        "snap_para:s"=>\$snap_para,
        "fgenesh"=>\$fgenesh,
        "geneid"=>\$geneid,
        "geneid_para:s"=>\$geneid_para,
        "prefix_for_Denovo:s"=>\$prefix_for_Denovo,
        "Dcpu:i"=>\$Dcpu,
        "Drun:s"=>\$Drun,
#        "shape:s"=>\$Shape,
    "EST"=>\$EST,
        "EST_species:s"=>\@EST_species,
        "input_EST:s"=>\@input_EST,
        "Ecpu:i"=>\$Ecpu,
        "Erun:s"=>\$Erun,
        "identity:s"=>\$Identity_cutoff,
        "alignrate:s"=>\$Alignrate_cutoff,
        "dbcut:i"=>\$DBcut,
        "cluster:s"=>\$Cluster,
        "tophit:i"=>\$Best_hit,
#       "prefix_for_EST:s"=>\$Prefix_for_EST,
    "cDNA"=>\$cDNA,
        "cDNA_species:s"=>\@cDNA_species,
        "input_cDNA:s"=>\@input_cDNA,
#       "prefix_for_cDNA:s"=>\$Prefix_for_cDNA,
    "Glean"=>\$Glean,
        "genome_name:s"=>\$genome_name,
        "parameter_yaml:s"=>\$parameter_yaml,
#       "lines_for_glean:i"=>\$lines_for_glean,
        "minlen:i"=>\$minlen,
        "minintron:i"=>\$minintron,
        "maxintron:i"=>\$maxintron,
        "cds:s"=>\$cds,
        "scaf:s"=>\$scaf,
        "resource_glean:s"=>\$Resource_glean,
        "cpu_for_glean:i"=>\$cpu_for_glean,
	"Process"=>\$Process,
		"SP_cuts:i"=>\$SP_cuts,
		"SP_cpu:i"=>\$SP_cpu,
		"blast_evalue:s"=>\$blast_evalue,
		"blast_vf:s"=>\$blast_vf,
		"rna_list:s"=>\$rna_list,
		"soap_vf:s"=>\$soap_vf,
		"Cover_cutoff:s"=>\$Cover_cutoff,
    "CCG"=>\$CCG,
        "reads_list:s"=>\$reads_list,
        "train_gff:s"=>\$train_gff,
        "combine_gff:s"=>\$gff_to_combine,
        "tsco:i"=>\$tscore,
        "tstart:i"=>\$tstart,
        "tstop:i"=>\$tstop,
	"qual_system:i"=>\$qual_system,
        "alignment:s"=>\$alignment,
        "build_vf:s"=>\$build_vf,
        "align_vf:s"=>\$align_vf,
        "max_i_len:i"=>\$max_i_len,
        "asmb_vf:s"=>\$asmb_vf,
        "asmb_proc_num:i"=>\$asmb_proc_num,
        "min_cds_len:i"=>\$min_cds_len,
        "orf:s"=>\$CCG_orf,
        "CCG_step:s"=>\$CCG_step,
        "CCG_name:s"=>\$CCG_name,
    "Verbose"=>\$Verbose,
    "Help"=>\$Help,
    "queue:s"=>\$Queue,
    "pro_code:s"=>\$Pro_code,
);

die `pod2text $0` if(@ARGV==0 || $Help);

my $config_file = "$Bin/../../config.txt";
my $input_fa = shift;
my $remask_genome = shift;
my $seq_name = basename($input_fa);
my $remask_name=basename($remask_genome);
my $dir = `pwd`; chomp($dir);
my $species_name=$1 if ($seq_name=~/^(\w+)\./);
$genome_name ||= "whole";

$opt.="--queue $Queue " if (defined $Queue);
$opt.="--pro_code $Pro_code " if (defined $Pro_code);

### Part 1 ### Homolog
if (defined $Homolog)
{
	$evidence .= "Homolog_";
	$Hcpu ||= 50;
	$Hrun ||= "qsub";
	$Blast_eval ||= "1e-5";
	$Align_rate ||= 0.25;
    $Extend_len ||= 500;
	$Step ||= "1234";
	$Lines ||= 500;
    
    my $protein_map_genome;
    if(defined $Genblasta){
        $Genblasta_opt ||="Genblasta_opt=\"-p T -e 1e-2 -g T -f F -a 0.5 -d 100000 -r 100 -c 0.5 -s -100\"";
        $Filter_rate ||= 0.7;
        $Hresource ||= "vf=1.5G";
        
        $protein_map_genome = parse_config($config_file,"protein_map_genome_genBlastA");
        $protein_map_genome.=" --resource $Hresource " if (defined $Hresource);
        $protein_map_genome.=" --genblasta_opt $Genblasta_opt " if (defined $Genblasta_opt);
        $protein_map_genome.=" --filter_rate $Filter_rate " if (defined $Filter_rate);
        $protein_map_genome.=" --net $Net " if (defined $Net);
        $protein_map_genome.=" --rgene $Ref_gene " if (defined $Ref_gene);
    }else{
        $protein_map_genome = parse_config($config_file,"protein_map_genome");
        $protein_map_genome.=" --tophit $Htophit " if (defined $Htophit);
        $protein_map_genome.=" --blast_eval $Blast_eval " if (defined $Blast_eval);
    }
    
	my $formatdb = parse_config($config_file,"formatdb");
	my @Inputpep_name = ();
	open(OUT , ">" . "STEP01-1_homolog_work.sh") or die $!;
	print OUT "mkdir Homolog\ncd Homolog\n";
    my $genome_seq;
    if(defined $Hremask){
        $genome_seq=$remask_genome;
    }else{
        $genome_seq=$input_fa;
    }
#    print OUT "$formatdb -i $genome_seq -p F -o T ;\n";
	for ($i = 0; $i < @Species; $i++)
	{
		print OUT "mkdir $Species[$i]\ncd $Species[$i]\n";
        print OUT "ln -s $genome_seq\n";
        my $basename=`basename $genome_seq`;
        chomp $basename;
		print OUT "nohup perl $protein_map_genome $opt --cpu $Hcpu --run $Hrun --align_rate $Align_rate --extend_len $Extend_len --step $Step --lines $Lines  --verbose $Inputpep[$i] $dir/Homolog/$Species[$i]/$basename &\n";
		print OUT "cd ..\n";
		$Inputpep_name[$i] = basename($Inputpep[$i]);
        if($protein_map_genome=~/genBlastA/){
            $all_glean_gff .= "--homolog $dir/Homolog/$Species[$i]/$Inputpep_name[$i].genblast.genewise.filter.gff ";
            $all_stat .= "Homolog/$Species[$i] $dir/Homolog/$Species[$i]/$Inputpep_name[$i].genblast.genewise.filter.gff\n";
        }else{
            $all_glean_gff .= "--homolog $dir/Homolog/$Species[$i]/$Inputpep_name[$i].solar.genewise.gff ";
            $all_stat .= "Homolog/$Species[$i] $dir/Homolog/$Species[$i]/$Inputpep_name[$i].solar.genewise.gff\n";
        }
	}
	print OUT "cd ..\n";
	close(OUT);
}

### Part 2 ### Denovo
if (defined $Denovo)
{
	my $denovo_gene_prediction = parse_config($config_file,"denovo_gene_prediction"); #run augustus, genescan, fgenesh
	my $GlimmerHMM = parse_config($config_file,"GlimmerHMM_prediction");
	my $Snap = parse_config($config_file,"snap_prediction");
	my $GeneID = parse_config($config_file,"geneid_prediction");
    my $fgenesh_para = parse_config($config_file,"fgenesh_para");
	$species_for_augustus ||= "arabidopsis";
    $genescan_para ||= parse_config($config_file,"genscan_para");
    $glimmerhmm_para ||= parse_config($config_file,"glimmerHMM_para");
    $snap_para ||= parse_config($config_file,"snap_para");
    $geneid_para ||= parse_config($config_file,"geneid_para");
	$Dcpu ||= 20;
	$Drun ||= "qsub";
#	$Shape ||= "linear";
	$Denovo_parameter = "";
	$Denovo_parameter .= " --prefix $prefix_for_Denovo " if(defined $prefix_for_Denovo);
	open(OUT , ">" . "STEP01-2_denovo_work.sh") or die $!;
	print OUT "mkdir Denovo\ncd Denovo\n";
	if (defined $augustus)
	{
		print OUT "mkdir Augustus\ncd Augustus\n";
		$evidence .= "Augustus_";
		print OUT "nohup perl $denovo_gene_prediction $opt --augustus $species_for_augustus --cpu $Dcpu $remask_genome &\n";
		print OUT "cd ..\n";
		$all_glean_gff .= "--denovo $dir/Denovo/Augustus/$remask_name.augustus.gff.check.gff ";
		$all_stat .= "Denovo/Augustus $dir/Denovo/Augustus/$remask_name.augustus.gff.check.gff\n";
	}
	if (defined $genescan)
	{
		print OUT "mkdir Genescan\ncd Genescan\n";
		$evidence .= "Genscan_";
		print OUT "nohup perl $denovo_gene_prediction $opt --genscan $genescan_para $Denovo_parameter --cpu $Dcpu --run $Drun $remask_genome &\n";
		print OUT "cd ..\n";
		$all_glean_gff .= "--denovo $dir/Denovo/Genescan/$remask_name.genscan.gff.check.gff ";
		$all_stat .= "Denovo/Genescan $dir/Denovo/Genescan/$remask_name.genscan.gff.check.gff\n";
	}

	if (defined $fgenesh)
	{
		print OUT "mkdir Fgenesh\ncd Fgenesh\n";
		$evidence .= "Fgenesh_";
		print OUT "nohup perl $denovo_gene_prediction $opt --fgenesh $fgenesh_para $Denovo_parameter --cpu $Dcpu --run $Drun $remask_genome &\n";
		print OUT "cd ..\n";
		$all_glean_gff .= "--denovo $dir/Denovo/Fgenesh/$remask_name.fgenesh.gff.check.gff ";
		$all_stat .= "Denovo/Fgenesh $dir/Denovo/Fgenesh/$remask_name.fgenesh.gff.check.gff\n";
	}

	if (defined $glimmerhmm)
	{
		print OUT "mkdir GlimmerHMM\ncd GlimmerHMM\n";
		$evidence .= "GlimmerHMM_";
		print OUT "nohup perl $GlimmerHMM $opt --glimmerHMM $glimmerhmm_para --cpu $Dcpu --run $Drun $remask_genome &\n";
		print OUT "cd ..\n";
		my @tmp_arr = ();
		@tmp_arr = split("/",$glimmerhmm_para);
		my $num = @tmp_arr - 1;
		$all_glean_gff .= "--denovo $dir/Denovo/GlimmerHMM/$remask_name.glimmerHMM.change.gff ";
		$all_stat .= "Denovo/GlimmerHMM $dir/Denovo/GlimmerHMM/$remask_name.glimmerHMM.change.gff\n";
	}
	
    if (defined $snap)
	{
		print OUT "mkdir SNAP\ncd SNAP\n";
		$evidence .= "SNAP_";
		print OUT "nohup perl $Snap $opt --snap $snap_para --cpu $Dcpu --run $Drun --cutf $Dcpu  $remask_genome &\n";
		print OUT "cd ..\n";
		$all_glean_gff .= "--denovo $dir/Denovo/SNAP/$remask_name.snap.gff ";
		$all_stat .= "Denovo/SNAP $dir/Denovo/SNAP/$remask_name.snap.gff\n";
	}

	if (defined $geneid)
	{
		print OUT "mkdir GeneID\ncd GeneID\n";
        $evidence .= "GeneID_";
		print OUT "nohup perl $GeneID $opt -geneid $geneid_para --cpu  $Dcpu --run $Drun --cutf $Dcpu $remask_genome &\n";
		print OUT "cd ..\n";
        $all_glean_gff .= "--denovo $dir/Denovo/GeneID/$remask_name.geneid.gff ";
        $all_stat .= "Denovo/GeneID $dir/Denovo/GeneID/$remask_name.geneid.gff\n";
	}
	print OUT "cd ..\n";
	close(OUT);
}

### Part 3 ### EST
if (defined $EST)
{
	$evidence .= "EST_";
	my $est_map_genome = parse_config($config_file,"est_map_genome");
	$Ecpu ||= 100;
	$Erun ||= "qsub";
	$Identity_cutoff ||= 0.95;
	$Alignrate_cutoff ||= 0.95;
	$DBcut ||= 3;
	$Cluster ||='T';
	open(OUT , ">" . "STEP01-3_EST_work.sh") or die $!;
	print OUT "mkdir EST\ncd EST\n";
	for ($i = 0; $i < @EST_species; $i++)
	{
		$EST_parameter = "";
		$EST_parameter .= " --tophit $Best_hit " if(defined $Best_hit);
		$EST_parameter .= " --prefix $EST_species[$i] ";
		print OUT "mkdir $EST_species[$i]\ncd $EST_species[$i]\n";
		print OUT "nohup perl $est_map_genome $opt --cpu $Ecpu --run $Erun --identity $Identity_cutoff --alignrate $Alignrate_cutoff --dbcut $DBcut --cluster $Cluster $EST_parameter $input_EST[$i] $input_fa &\n";
		print OUT "cd ..\n";
		$input_EST_name[$i]=basename($input_EST[$i]);
		$all_glean_gff .= "--EST $dir/EST/$EST_species[$i]/$input_EST_name[$i].blat.sim4.gff.pasa.gff ";
		$all_stat .= "EST/$EST_species[$i] $dir/EST/$EST_species[$i]/$input_EST_name[$i].blat.sim4.gff.pasa.gff\n";
	}
	print OUT "cd ..\n";
	close(OUT);
}
	 ### cDNA
if (defined $cDNA)
{
	$evidence .= "EST_";
	my $cDNA_map_genome = parse_config($config_file,"cDNA_map_genome");
	$Ecpu ||= 100;
	$Erun ||= "qsub";
	$Identity_cutoff ||= 0.95;
	$Alignrate_cutoff ||= 0.95;
	$DBcut ||= 10;
	open(OUT , ">" . "STEP01-4_cDNA_work.sh") or die $!;
	print OUT "mkdir cDNA\ncd cDNA\n";
	for ($i = 0; $i < @cDNA_species; $i++)
	{
		my $cDNA_parameter = "";
		$cDNA_parameter .= " --tophit $Best_hit " if(defined $Best_hit);
		$cDNA_parameter .= " --prefix $cDNA_species[$i] ";
		print OUT "mkdir $cDNA_species[$i]\ncd $cDNA_species[$i]\n";
		print OUT "nohup perl $cDNA_map_genome $opt --cpu $Ecpu --run $Erun --identity $Identity_cutoff --alignrate $Alignrate_cutoff --dbcut $DBcut $cDNA_parameter $input_cDNA[$i] $input_fa &\n";
		print OUT "cd ..\n";
		$input_cDNA_name[$i]=basename($input_cDNA[$i]);
		$all_glean_gff .= "--cDNA $dir/cDNA/$cDNA_species[$i]/$input_cDNA_name[$i].blat.sim4.gff ";
		$all_stat .= "cDNA/$cDNA_species[$i] $dir/cDNA/$cDNA_species[$i]/$input_cDNA_name[$i].blat.sim4.gff\n";
	}
	print OUT "cd ..\n";
	close(OUT);
}

### Part 4 ### GLEAN
if (defined $Glean)
{
	my $glean_program = parse_config($config_file,"glean_program");
	my $check_parameter_yaml = parse_config($config_file,"check_parameter_yaml");
	$parameter_yaml ||= parse_config($config_file,"parameter_yaml");
	$run_for_glean ||= "qsub";
	#$lines_for_glean ||= 100;
	$minlen ||= 150;
	$scaf ||= 0;
	$cds ||= 0;
	$minintron ||= 11;
	$maxintron ||= 10000;
	$cpu_for_glean ||= 30;
	$Resource_glean ||= "vf=0.3G";

	open(OUT , ">" . "STEP02-1_glean_work.sh") or die $!;
	print OUT "mkdir glean\ncd glean\n";
	unless (defined $evidence) {
		`perl $Bin/File_find.pl $dir >search`;
		open IN,'search' or die $!;
		$evidence=<IN>;
		chomp $evidence;
		$all_glean_gff=<IN>;
		chomp $all_glean_gff;
		close IN;
	}
	#print OUT "cp $glean_par_file $dir/glean/parameter.yaml\n";
	print OUT "perl $check_parameter_yaml $evidence $parameter_yaml >$dir/glean/parameter.yaml \n";
	print OUT "nohup perl $glean_program $opt --genome_name $genome_name --genome $input_fa -minlen $minlen -minintron $minintron -maxintron $maxintron -run $run_for_glean -cpu $cpu_for_glean --YAML $dir/glean/parameter.yaml --cds $cds --scaf $scaf --resource $Resource_glean $all_glean_gff &\n";
	print OUT "cd ..\n";
	$all_stat .= "Glean $dir/glean/result/$species_name.gene.gff\n";
	close(OUT);
	
}
### Part 5 ### Process
if (defined $Process)
{
	my $auto_glean_process = parse_config($config_file,"auto_glean_process");
	$SP_cuts ||= 100;
	$SP_cpu ||= 30;
	$blast_evalue ||= 1e-5;
	$blast_vf ||= "vf=1G";
	$soap_vf ||= "vf=4G";
	$Cover_cutoff ||= 80;

	open(OUT , ">" . "STEP02-2_glean_process.sh") or die $!;
	print OUT "mkdir process\ncd process\n";
	my @all_glean = split (/--/, $all_glean_gff);
	my $homo_gff;
	foreach my $i ( 1 .. $#all_glean){
		my @tmp = split (/\s+/, $all_glean[$i]);
		$homo_gff .= "--homo_gff $tmp[1] " if ($tmp[0] eq "homolog");
	}
	print OUT "nohup perl $auto_glean_process $opt --homolog $homo_gff --swissprot --cuts $SP_cuts --cpu $SP_cpu --evalue $blast_evalue --blast_vf $blast_vf --RNAseq --rna_list $rna_list --soap_vf $soap_vf --cutoff $Cover_cutoff $input_fa $dir/glean/result/$species_name.gene.gff &\n";
	$all_stat .= "Process $dir/process/$species_name.gene.process.gff\n";
	print OUT "cd ..\n";
}

	 ### CCG
if (defined $CCG){
	my $CCG_program = parse_config($config_file,"auto_CCG");
	my $training_markov5 = parse_config($config_file,"training_markov5");
	my $perfect_gene= parse_config($config_file,"perfect_gene");
	
	unless (defined $train_gff){
#		$all_glean_gff=~/--homolog\s+(\S+)/ || $all_glean_gff=~/--\S+\s+(\S+)/;
        my $all_gff_evidence=$all_glean_gff;
        $all_gff_evidence=~s/^--//g;
		my @all_gff_evidence=split /--/,$all_gff_evidence;
        foreach my $temp(@all_gff_evidence){
             if($temp=~/(\S+)\s+(\S+)/ && $1 eq "homolog"){
#			 if($temp=~/(\S+)\s+(\S+)/ && $1 ne "denovo"){
                 $train_gff.=" $2 ";
             }
        }
	}

    $tscore ||="95";
    $tstart ||="10";
    $tstop ||="10";
    	$qual_system ||= 64;
	$alignment ||="hisat";
	$build_vf ||="vf=4G";
	$align_vf ||="vf=5G";
	$max_i_len ||=500000;
	$asmb_vf ||="vf=5G";
	$asmb_proc_num ||=1;
	$min_cds_len ||=150;
	$CCG_orf ||="T";
	$CCG_step ||="123456";
	$CCG_name ||="CCG";
	$reads_list ||= $rna_list;
	$gff_to_combine ||="$dir/process/$species_name.gene.process.gff";

	open OUT,">STEP02-3_CCG_work.sh" || die $!;
	print OUT "mkdir CCG\ncd CCG\n";
	print OUT "mkdir Markov5\ncd Markov5\n";
#	my $genome_for_perfect_gene=(defined $remask_genome)? $remask_genome : $input_fa;
	print OUT "perl $perfect_gene --sco $tscore --start $tstart --stop $tstop $remask_genome $train_gff\n";
	my $pefect_result=basename $remask_genome;
	$pefect_result="$pefect_result.gff.nr.gff";
	## To choose 1000 perfect genes??
	print OUT "# the work of training markov5 need to be qsub ...\n";   # add on 20161107
	print OUT "perl $training_markov5 $pefect_result $remask_genome > training_markov5.log\n";
	my $markov5_para = "$dir/CCG/Markov5/markov_5.param";
	print OUT "cd ..\n";
    print OUT "mkdir Data\ncd Data\n";
    my $CCG_input_fa;
    if ($input_fa=~/\.fa$/){
        print OUT "ln -s $input_fa\n";
        my $new_name=`basename $input_fa`;
        chomp $new_name;
        $CCG_input_fa="$dir/CCG/Data/$new_name";
    }elsif ($input_fa=~/\.fasta$/){
        my $new_name=`basename $input_fa`;
        chomp $new_name;
        $new_name=~s/fasta/fa/g;
        print OUT "ln -s $input_fa $new_name\n";
        $CCG_input_fa="$dir/CCG/Data/$new_name";
    }else{
        my $new_name=`basename $input_fa`;
        chomp $new_name;
        $new_name=$new_name.".fa";
        print OUT "ln -s $input_fa $new_name\n";
        $CCG_input_fa="$dir/CCG/Data/$new_name";
    }
    print OUT "cd ..\n";
	print OUT "nohup perl $CCG_program $opt -qual_system $qual_system -alignment hisat -build_vf $build_vf -align_vf $align_vf -max_i_len $max_i_len -asmb_vf $asmb_vf -asmb_proc_num $asmb_proc_num -min_cds_len $min_cds_len -orf $CCG_orf -step $CCG_step -species $CCG_name $CCG_input_fa $reads_list $markov5_para $gff_to_combine &\n"; ####
	$all_stat .= "CCG $dir/CCG/CCG.rebuild.gff.reorder.noisoforms.gff\n";
}

open (OUT ,">STEP03_gene_statistics.sh") or die $!;
my $gene_set_statistics="$genome_name\.gene\.evidence_statistics\.xls";
print OUT "mkdir statistics\ncd statistics\n";
$all_stat=~s/\n*$//;
my @gene_sets=split /\n/,$all_stat;
my $stat_gene=parse_config($config_file,"stat_gene");
my $sta_tag=0;
foreach my $gene_set (@gene_sets) {
	my @sin_set=split /\s/,$gene_set;
	print OUT "ln -s $sin_set[1]\n";
	if ($sta_tag==0){
		print OUT "perl $stat_gene -name $sin_set[0] $sin_set[1] >>$gene_set_statistics\n";
		$sta_tag++;
	}else{
		print OUT "perl $stat_gene -name $sin_set[0] -notag $sin_set[1] >>$gene_set_statistics\n";
	}
}
close OUT;
open(OUT , ">" . "STEP04_delete_tmp_files.sh") or die $!;
if (defined $Homolog) {
	foreach my $Specie (@Species) {
		print OUT "mv $dir/Homolog/$Specie/*gff $dir/Homolog/\n";
		print OUT "rm -r ./Homolog/$Specie/*\n";
		print OUT "mv $dir/Homolog/*gff $dir/Homolog/$Specie/\n";
	}
}
if (defined $augustus) {
	print OUT "mv $dir/Denovo/Augustus/*gff $dir/Denovo/\n";
	print OUT "rm -r $dir/Denovo/Augustus/*\n";
	print OUT "mv $dir/Denovo/*gff $dir/Denovo/Augustus/\n";
}
if (defined $genescan) {
	print OUT "mv $dir/Denovo/Genescan/*gff $dir/Denovo/\n";
	print OUT "rm -r $dir/Denovo/Genescan/*\n";
	print OUT "mv $dir/Denovo/*gff $dir/Denovo/Genescan/\n";
}
if (defined $glimmerhmm) {
	print OUT "mv $dir/Denovo/GlimmerHMM/*gff $dir/Denovo/\n";
	print OUT "rm -r $dir/Denovo/GlimmerHMM/*\n";
	print OUT "mv $dir/Denovo/*gff $dir/Denovo/GlimmerHMM/\n";
}
if (defined $snap) {
	print OUT "mv $dir/Denovo/SNAP/*gff $dir/Denovo/\n";
	print OUT "rm -r $dir/Denovo/SNAP/*\n";
	print OUT "mv $dir/Denovo/*gff $dir/Denovo/SNAP/\n";
}
if (defined $EST) {
	foreach my $EST_specie (@EST_species) {
		print OUT "mv $dir/EST/$EST_specie/*pasa.gff $dir/EST/\n";
		print OUT "rm -r $dir/EST/$EST_specie/*\n";
		print OUT "mv $dir/EST/*pasa.gff $dir/EST/$EST_specie/\n";
	}
}
if (defined $cDNA) {
	foreach my $cDNA_specie (@cDNA_species) {
		print OUT "mv $dir/cDNA/$cDNA_specie/*gff $dir/cDNA/\n";
		print OUT "rm -r $dir/cDNA/$cDNA_specie/*\n";
		print OUT "mv $dir/cDNA/*gff $dir/cDNA/$cDNA_specie/\n";
	}
}
if (defined $Glean) {
	print OUT "rm -r $dir/glean/result/glean* $dir/glean/result/parameter* $dir/glean/result/lca/\n";
}
if (defined $CCG){
	print OUT "rm -r $dir/CCG/compare.sh.*.qsub $dir/CCG/stdout.tmap.split $dir/CCG/tophat.sh.14646.qsub/ \n";
}
close OUT;
