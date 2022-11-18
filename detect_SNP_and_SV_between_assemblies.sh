
# align the contigs of two assemblies and identify SNPs with MUMmer
/usr/bin/MUMmer3.23/nucmer -p nucmer --mum GH06.contigs ZY821.contigs
/usr/bin/MUMmer3.23/delta-filter -1 nucmer.delta  > nucmer.delta.filter
/usr/bin/MUMmer3.23/show-coords -rcl nucmer.delta.filter > nucmer.delta.filter.coords
/usr/bin/MUMmer3.23/show-snps -Clr -x 1 -T nucmer.delta.filter > nucmer.delta.filter.snps
python /usr/bin/PythonNGSTools-master/MUMmerSNPs2VCF.py nucmer.delta.filter.snps nucmer.delta.filter.snps.vcf

#Detect SV with Assemblytics
export PATH="/usr/bin/Assemblytics-master/:$PATH"
alias R='/usr/bin/R'
alias Rscript='/usr/bin/Rscript'
/usr/bin/Assemblytics-master/Assemblytics nucmer.delta.filter ./Bnapus 10000 /usr/bin/Assemblytics-master/
/usr/bin/SURVIVOR-master/Debug/SURVIVOR convertAssemblytics Bnapus.Assemblytics_structural_variants.bed 50 Bnapus.Assemblytics_structural_variants.min50bp.vcf

