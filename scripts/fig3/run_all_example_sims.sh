#!/bin/bash
set -e

bash ./fig3.sh smaller_pop_size/fig3.toml
bash ./fig3.sh smaller_pop_size/fig3_higher_vax_direct_effect.toml
bash ./fig3.sh smaller_pop_size/fig3_lower_vax_direct_effect.toml