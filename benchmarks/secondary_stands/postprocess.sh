#!/bin/bash

# See also the comment in benchmarks/secondary_stands/guess.ins about this test,
# about using it to test savestate/restart, and about other output of this benchmark.
#
# When the benchmark is run as part of a stardard run of the full benchmark suite:
# - it is run without restart.
# - the postprocess script (not implemented yet) checks that columns Natural and mixed 
#   in the file cmass_sts.out are identical up through 1950.
# - a fail to that test indicates that developer forgot to initiate something when stands are 
#   created later than year zero.
# - this is a works-or-breaks type of test, similarly to crop_mixed_sites - its html-report 
#   page will only state whether it passed or not.
# 
# When the benchmark is run as part of a test of the model's savestate/restart:
# - use the stateyear in the file benchmarks\secondary_stands\config\guess.ins.
# - a fail indicates, for example, that the variable-lists in serialize() functions are not 
#   complete, etc; or could indicate an error to cloning or running multiple stands. 


cat cmass_sts.out | tail -n+2 | head -n50\
 | awk '$4 != $5 {print "\nTEST FAILED! Benchmark test secondary_stands failed: Natural and mixed columns differ at line", NR, "in cmass_sts.out","\nThis indicates a problem with secondary stands (e.g. a member variable was not inititated in constructor, etc.)"; exit;}'\
 >testresult.txt

cat testresult.txt >>guess.log


# Benchmark report

describe_benchmark "LPJ-GUESS - Secondary stands technical test"

if [ -z "$(cat testresult.txt)" ]; then		# Test passed.
  echo "Test passed." >testresult.txt
fi						# If test failed, use the pre-existing text in testresult.txt

describe_textfile testresult.txt

