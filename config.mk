# Configurations for the TATSP runner

# Max solver times
MAXTIME_LB_LP = 120
MAXTIME_LB_RLXLAG = 30
MAXTIME_LB_COLGEN = 1030
MAXTIME_UB_LP = 60
MAXTIME_UB_RLXLAG = 60
MAXTIME_UB_COLGEN = 1060
MAXTIME_ILP = 120
# MAXTIME_ILP = 60

# Other parameters
SEED_NUMBER = 1234
RA_PARAM = 291111
# VERBOSE_MODE = --verbose

# Choose between directories to change instances
# INSTANCES_DIR = instances/instances_release_1
INSTANCES_DIR = instances/instances_release_2
# INSTANCES_DIR = instances

# Toy instance for quick testing
TEST_INSTANCE := instances/inst-slide8.txt
# TEST_INSTANCE := instances/instances_release_1/grf1.txt
# TEST_INSTANCE := instances/instances_release_2/grf101.txt
