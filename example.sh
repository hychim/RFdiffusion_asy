./scripts/run_inference_asy.py \
	diffuser.T=50 \
    inference.symmetry="C4" \
	inference.asy_motif="True" \
    inference.num_designs=1 \
    inference.output_prefix=./samples/auto/1DBB_0-0-180_25_5 \
	inference.input_pdb=./examples/input_pdbs/1DBB_interface_only.pdb \
    "contigmap.contigs=[5/Y471-485/15-35/X95-101/15-35/X246-251/15-35/X265-271/1/X273-278/15-35/X313-322/5/0]" \
	inference.motif_drag="True" \
	inference.asy_motif_weight=1 \
	inference.asy_motif_rot_range=[180,90,180] \
	inference.asy_motif_dist=25 \
	inference.asy_motif_dist_range=5 \
	inference.random_drag="True"
