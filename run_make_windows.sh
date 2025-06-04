
# minimal arguments with gtf
echo "---minimal arguments with gtf"
./make_windows.R --gtf data/example.gtf --upstream 1000 --downstream 2000 --out out_tolga/minimal_gtf.bed

# minimal arguments with bed
echo "---minimal arguments with bed"
./make_windows.R --bed data/example_annotation_1000up_2000down.bed --upstream 1000 --downstream 2000 --out out_tolga/minimal_bed.bed

# filter wrt feature="gene" biotype="lncRNA"
echo "---filter wrt feature=\"gene\" biotype=\"lncRNA\""
./make_windows.R --gtf data/example.gtf --upstream 1000 --downstream 2000 --out out_tolga/example_gtf_gene_lncRNA.bed --feature gene --biotype lncRNA

# edge case with edge make_windows
echo "---edge case with edge make_windows"
./make_windows.R --gtf data/example_edges.gtf --upstream 1000 --downstream 2000 --out out_tolga/example_gtf_gene_lncRNA_edges.bed --feature gene --biotype lncRNA

