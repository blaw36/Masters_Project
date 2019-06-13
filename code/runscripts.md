## WaveQTL_HMT running calls:

- All 24 SNPs on all 1024 phenotypes, using the pre-processing steps as per the WaveQTL software manual.
- Will use its "use/do not use" vector for WCs
- Will use quantile-transformed for likelihood ratio testing
- Will use non-quantile transformed for measuring effect size in data space.

### Tying
Several tying techniques:

#### Group 1: Top (root) of tree on its own level

1) Tree-level tying (uses default group file, which is based on tree levels)
```
../../WaveQTL -gmode 1 -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o tree_tie -f 1024 -hmt 1 > stdout.txt 2> stderr.txt 
```

2) Tree levels 2-3 tied together
```
../../WaveQTL -gmode 1 -group g23_1024.txt -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o g23_tie -f 1024 -hmt 1 > stdout.txt 2> stderr.txt 
```

3) Tree levels 2-4 tied together
```
../../WaveQTL -gmode 1 -group g24_1024.txt -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o g24_tie -f 1024 -hmt 1 > stdout.txt 2> stderr.txt 
```

4) Tree levels 2-5 tied together
```
../../WaveQTL -gmode 1 -group g25_1024.txt -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o g25_tie -f 1024 -hmt 1 > stdout.txt 2> stderr.txt 
```

#### Group 2: Top (root) of tree grouped with the first tying level

5) Tree levels 1-3 tied together
```
../../WaveQTL -gmode 1 -group g13_1024.txt -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o g13_tie -f 1024 -hmt 1 > stdout.txt 2> stderr.txt 
```

6) Tree levels 1-4 tied together
```
../../WaveQTL -gmode 1 -group g14_1024.txt -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o g14_tie -f 1024 -hmt 1 > stdout.txt 2> stderr.txt 
```

7) Tree levels 1-5 tied together
```
../../WaveQTL -gmode 1 -group g15_1024.txt -g ../../data/dsQTL/chr17.10160989.10162012.2kb.cis.geno -p WCs.txt -u use.txt -o g15_tie -f 1024 -hmt 1 > stdout.txt 2> stderr.txt 
```
