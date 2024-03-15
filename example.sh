./scripts/run_inference_asy.py \
	diffuser.T=50 \
        inference.symmetry="C4" \
	inference.asy_motif="True" \
        inference.num_designs=5 \
        inference.output_prefix=./output/test_z_ran10/z_ran10 \
	inference.input_pdb=./examples/input_pdbs/1DBB_interface_only.pdb \
        "contigmap.contigs=[5/Y471-485/20/X95-101/20/X246-251/20/X265-271/1/X273-278/20/X313-322/5/0]"\
