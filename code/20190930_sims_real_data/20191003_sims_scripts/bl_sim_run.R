source("/data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v4/bl_sim_functions.R")

# Test run for 1st dataset
effect_to_prop_and_sim(shrink_factor = 1
                       ,datasetNum = 1
                       ,effect.size.path = "/data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/data/smooth.ratio.3."
                       ,geno.path = "/data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/data/geno70.dat"
                       ,raw.dat.path = "/data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v3/data/pheno.dat"
                       ,output.data.path = "/data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v4/data/"
                       ,read.depth.ratio = 1
                       ,over.dispersion = 1/70/70/10
                       ,multipleSig = 1
                       ,wavelet.preprocess = TRUE
                       ,DESeq.preprocess = TRUE
                       ,wd.path = "/data/cephfs/punim0614/shared/shared_data/internal/multi.scale/multiseq/forBrendan/simulation_manyQTLfinal_v4/"
                       ,output.dir.name = "full")