./scripts/run_inference_asy.py \
	diffuser.T=50 \
        inference.symmetry="C4" \
	inference.asy_motif="True" \
        inference.num_designs=3 \
        inference.output_prefix=./output/test_ran20_randrag/test_ran20_randrag \
	inference.input_pdb=./examples/input_pdbs/1DBB_interface_only.pdb \
        "contigmap.contigs=[5/Y471-485/15-35/X95-101/15-35/X246-251/15-35/X265-271/1/X273-278/15-35/X313-322/5/0]"\
