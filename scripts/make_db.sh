#!/usr/bin/env bash
#SBATCH --cpus=4
#SBATCH --mem=10
. ~soft_cvi_114/initializes/init_fln
make_user_db.rb -u vertebrates -l -t Actinopterygii -n Actinopterygii_db
make_user_db.rb -u vertebrates -l -f -t Danio -n Danio_db
