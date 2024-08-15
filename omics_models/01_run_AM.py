import argparse
from alphamissense.data import pipeline_missense
from alphamissense.model import config
from alphamissense.model import modules_missense
import jax
import numpy as np
import json
import haiku as hk


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("sequenceFile", help="Sequence to be predicted in fasta format", type = str)
    parser.add_argument("proteinID", help="Name of target protein", type = str)
    parser.add_argument("referenceAA", help="Reference AA", type = str)
    parser.add_argument("position", help="AA position", type = int)
    parser.add_argument("targetAA", help="Target AA", type = str)
    # subfolder param where you want to place the embedding output
    # parser.add_argument("function", help="Mutational Function", type = str, default="gof")
    parser.add_argument("-db","--db_dir", help="directory to local MSA databases", type = str, default="DEFAULT PATH HERE")
    parser.add_argument("-o","--output", help="Output DIR", type = str)

    args = parser.parse_args()

    protein_sequence_file = args.sequenceFile
    DATABASES_DIR = args.db_dir

    with open('/path/to/protein_sequences.json', 'r') as file:
        gene_entry_dict = json.load(file)

    if args.proteinID in gene_entry_dict.keys():
        pid = gene_entry_dict[args.proteinID]
    else:
        pid = args.proteinID

    print(f"Embedding Gene: {args.proteinID}; Uniprot ID: {pid}")
    v_str = args.proteinID + f'_{args.referenceAA + str(args.position) + args.targetAA}'

    msa_dir = f'/output/path/to/msa_{v_str}'

    pipeline = pipeline_missense.DataPipeline(
        jackhmmer_binary_path="/path/to/jackhammer/bin/jackhmmer",  # Typically '/usr/bin/jackhmmer'.
        protein_sequence_file=protein_sequence_file,
        uniref90_database_path=DATABASES_DIR + '/uniref90/uniref90.fasta',
        mgnify_database_path=DATABASES_DIR + '/mgnify/mgy_clusters_2022_05.fa',
        small_bfd_database_path=DATABASES_DIR + '/small_bfd/bfd-first_non_consensus_sequences.fasta',
    )

    sample = pipeline.process(
        protein_id=pid,  # Sequence identifier in the FASTA file.
        reference_aa=args.referenceAA,  # Single capital letter, e.g. 'A'.
        alternate_aa=args.targetAA,
        position=args.position,  # Integer, note that the position is 1-based!
        msa_output_dir=msa_dir,
    )

    def _forward_fn(batch):
        model = modules_missense.AlphaMissense(config.model_config().model)
        return model(batch, is_training=False, return_representations=True)

    random_seed = 0
    prng = jax.random.PRNGKey(random_seed)
    params = hk.transform(_forward_fn).init(prng, sample)
    apply = jax.jit(hk.transform(_forward_fn).apply)
    output = apply(params, prng, sample)
    repr_out = f"./embs/embedding_{v_str}.npy"

    with open(repr_out, 'wb') as g:
        print(f"Saving pair representation embedding to {repr_out}")
        pair_rep = np.asarray(output["representations"]["pair"])
        np.save(g, pair_rep)