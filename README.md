collapse_to_uniq_seq.pl
=======================

Usage: /home/tristan/bin/collapse_to_uniq_seq.pl [options] <alignment> <seq output> <convertion table output>

Will collapse identical sequences. Ns and leading and ending gaps are not considered as differences. BioPerl is required.

Options:
        -h
                print this help
        -format: fasta|phylip|...
                give output format (default=fasta)
        -noamb, will convert all the ambiguity code into missing data (only for DNA data)
