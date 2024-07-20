## Adapted from utils/inference_utils.py
import numpy as np
import torch
import prody as pr
import sys
import os

import esm
from esm import FastaBatchedDataset, pretrained
from datasets.parse_chi import aa_idx2aa_short, get_onehot_sequence

args = sys.argv

protein_path = args[1]
protein_name = os.path.splitext(os.path.basename(protein_path))[0]

def get_sequences_from_pdbfile(file_path):
    sequence = None

    pdb = pr.parsePDB(file_path)
    seq = pdb.ca.getSequence()
    one_hot = get_onehot_sequence(seq)

    chain_ids = np.zeros(len(one_hot))
    res_chain_ids = pdb.ca.getChids()
    res_seg_ids = pdb.ca.getSegnames()
    res_chain_ids = np.asarray([s + c for s, c in zip(res_seg_ids, res_chain_ids)])
    ids = np.unique(res_chain_ids)

    for i, id in enumerate(ids):
        chain_ids[res_chain_ids == id] = i

        s_temp = np.argmax(one_hot[res_chain_ids == id], axis=1)
        s = ''.join([aa_idx2aa_short[aa_idx] for aa_idx in s_temp])

        if sequence is None:
            sequence = s
        else:
            sequence += (":" + s)

    return [sequence]

def compute_ESM_embeddings(model, alphabet, labels, sequences):
    # settings used
    toks_per_batch = 4096
    repr_layers = [33]
    include = "per_tok"
    truncation_seq_length = 1022

    dataset = FastaBatchedDataset(labels, sequences)
    batches = dataset.get_batch_indices(toks_per_batch, extra_toks_per_seq=1)
    data_loader = torch.utils.data.DataLoader(
        dataset, collate_fn=alphabet.get_batch_converter(truncation_seq_length), batch_sampler=batches
    )

    assert all(-(model.num_layers + 1) <= i <= model.num_layers for i in repr_layers)
    repr_layers = [(i + model.num_layers + 1) % (model.num_layers + 1) for i in repr_layers]
    embeddings = {}

    with torch.no_grad():
        for batch_idx, (labels, strs, toks) in enumerate(data_loader):
            #print(f"Processing {batch_idx + 1} of {len(batches)} batches ({toks.size(0)} sequences)")
            if torch.cuda.is_available():
                toks = toks.to(device="cuda", non_blocking=True)

            out = model(toks, repr_layers=repr_layers, return_contacts=False)
            representations = {layer: t.to(device="cpu") for layer, t in out["representations"].items()}

            for i, label in enumerate(labels):
                truncate_len = min(truncation_seq_length, len(strs[i]))
                embeddings[label] = representations[33][i, 1: truncate_len + 1].clone()
    return embeddings
    
print("Generating ESM language model embeddings..")

model_location = "esm2_t33_650M_UR50D"
model, alphabet = pretrained.load_model_and_alphabet(model_location)
model.eval()

if torch.cuda.is_available():
	model = model.cuda()

protein_sequence = get_sequences_from_pdbfile(protein_path)
labels, sequences = [], []

s = protein_sequence[0].split(':')
sequences.extend(s)

labels.extend([protein_name + '_chain_' + str(j) for j in range(len(s))])

lm_embeddings = compute_ESM_embeddings(model, alphabet, labels, sequences)

output_lm_embeddings = [lm_embeddings[f'{protein_name}_chain_{j}'] for j in range(len(s))]

torch.save(output_lm_embeddings, f"data/protein_embeddings/{protein_name}.pt")