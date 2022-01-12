"""Some basic testing and demo code for the bito module.

This script demonstrates the flag options for passing arguments to gradient and likelihood functions.

If you want to see the results of the print statements, use `pytest -s`.
"""

import sys
import json
import pprint
import pytest
import numpy as np
import bito
import bito.beagle_flags as beagle_flags
import sys

def gradients_with_flags_demo():
    # build instance
    inst = bito.rooted_instance("cheese")
    inst.read_newick_file("data/fluA.tree")
    inst.read_fasta_file("data/fluA.fa")
    inst.parse_dates_from_taxon_names(True)
    spec = bito.PhyloModelSpecification(
        substitution="GTR", site="weibull+4", clock="strict"
    )
    inst.prepare_for_phylo_likelihood(spec, 1, [beagle_flags.VECTOR_SSE], False)

    # these are flags/options for gradient and likelihood functions.
    import bito.phylo_flags as flags

    # these are keys that correspond to parts of the phylo model block map.
    import bito.phylo_keys as keys

    # NOTE: flags and keys that have the same name have the same value.
    # e.g. keys.SUBSTITUTION_MODEL == flags.SUBSTITUTION_MODEL

    # initialize model parameters (with enums)
    phylo_model_param_block_map = inst.get_phylo_model_param_block_map()
    phylo_model_param_block_map[keys.SUBSTITUTION_MODEL_RATES][:] = np.repeat(1 / 6, 6)
    phylo_model_param_block_map[keys.SUBSTITUTION_MODEL_FREQUENCIES][:] = np.repeat(
        1 / 4, 4
    )
    phylo_model_param_block_map[keys.SITE_MODEL_PARAMETERS][:] = np.array([0.5])
    phylo_model_param_block_map[keys.CLOCK_MODEL_RATES][:] = np.array([0.001])

    # Request and calculate gradients for RATIOS_ROOT_HEIGHT,
    # SUBSTITUTION_MODEL_FREQUENCIES, SUBSTITUTION_MODEL_RATES
    bito_grad = inst.phylo_gradients(
        # explicit flags: For flags that don't have associated values, can just be a list.
        [
            flags.RATIOS_ROOT_HEIGHT,
            flags.SUBSTITUTION_MODEL_FREQUENCIES,
            flags.SUBSTITUTION_MODEL_RATES,
        ],
        # run_with_default_flags:
        # (1) If set to true, all fields of phylo_model_block_map are populated unless overriden by an explicit flag.
        # (2) If set to false, no fields of phylo_model_block_map are populated unless overriden by an explicit flag.
        False,
    )[0]

    print(f"bito_grad_keys: {bito_grad.gradient.keys()}")
    gtr_rates_grad = np.array(bito_grad.gradient[keys.SUBSTITUTION_MODEL_RATES])
    gtr_freqs_grad = np.array(bito_grad.gradient[keys.SUBSTITUTION_MODEL_FREQUENCIES])
    ratios_root_height = np.array(bito_grad.gradient[keys.RATIOS_ROOT_HEIGHT])

    print(f"GTR rates gradient: {gtr_rates_grad}\n")
    print(f"GTR freqs gradient: {gtr_freqs_grad}\n")
    print(f"root height gradient: {ratios_root_height}\n")

    # This should trigger an exception because CLOCK_MODEL_RATES was not flagged to be computed, so the key does not exist in the output map.
    try:
        clock_grad = np.array(bito_grad.gradient[keys.CLOCK_MODEL_RATES])
        print("CHECK_THROW Failed: Key error not caught.")
    except:
        print("CHECK_THROW Successful: Error successfully caught.")
        print("ERROR: ", sys.exc_info()[0], "occurred.")

    # Above works if only boolean flags are used, otherwise use ordered tuples:
    bito_grad = inst.phylo_gradients(
        # explicit flags: For SET flags that require value, use ordered tuples. Non-SET flags just take boolean.
        [(flags.SET_GRADIENT_DELTA, 5.0)],
        # run_with_default_flags
        True,
    )[0]

    gtr_rates_grad = np.array(bito_grad.gradient[keys.SUBSTITUTION_MODEL_RATES])
    gtr_freqs_grad = np.array(bito_grad.gradient[keys.SUBSTITUTION_MODEL_FREQUENCIES])
    ratios_root_height = np.array(bito_grad.gradient[keys.RATIOS_ROOT_HEIGHT])

    print("GTR rates gradient: \n", gtr_rates_grad)
    print("GTR freqs gradient: \n", gtr_freqs_grad)
    print('root height gradient: \n', ratios_root_height)

    return inst

# run tests if called directly
if __name__ == "__main__":
    gradients_with_flags_demo()
