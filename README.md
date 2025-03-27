# Viral Syndromic Target Enrichment Pipeline (ViralSTEP)

Files and methodology pertaining to the development of virus hybridization probes or baits. VHost-Classifier and BaitsTools (v1.6.2) provided as submodules.

Help screen below **(Note: functions with appended "(WIP)" are still in development)**:

```
usage: viralstep <command> [<args>]

Accepted commands (add -h/--help after command for more details):
        ---------------------------------------------------------------------------------------
        Download Software, Packages, Public Databases
        ---------------------------------------------------------------------------------------
        download         Download all required software, lanuage interpreters, and packages (i.e., R, Python, Ruby, Perl) (WIP)
        downloadvhost    Download VHost-Classifier and update corresponding Virus-Host DB file (WIP)
        downloadbaits    Download required software and packages for bait tiling and filtering (i.e., BaitsTools) (WIP)
        downloadblastdb  Download (or update existing) BLASTn database (only nt)
        downloadstats    Download required software and packages for summary statistics and plots (WIP)

        ---------------------------------------------------------------------------------------
        Classifcation
        ---------------------------------------------------------------------------------------
        updatetax        Update required taxonomy files for Taxonkit (will install to $HOME/.taxonkit)
        generatedb       From a given FASTA file, generate classification database file (.csv)
        classify         Using an existing database file (.csv), create a separate database file (.csv) for a given FASTA file

        ---------------------------------------------------------------------------------------
        VHost-Classifier (https://github.com/Kzra/VHost-Classifier; PMID: 30824917)
        ---------------------------------------------------------------------------------------
        vhost_update     Update Virus-Host DB file (will default to VHost-Classifier location)
        vhost_filter     Invoke VHost-Classifier for a given FASTA file (filters for Mammalian)

        ---------------------------------------------------------------------------------------
        Bait Design
        ---------------------------------------------------------------------------------------
        baits            Run default bait design protocol for a given FASTA file (SKIPS syndromicbaits and humanwgc) (WIP)
        tilebaits        Run BaitsTools (tilebaits) for given input FASTA, bait length, offset, GC Content range
        tmbaits          Run melting temperature filter for given FASTA and temperature range (inclusive)
        syndromicbaits   Run syndromic filter for given FASTA, classification database file, and desired syndrome
        crosshybaits     Run cross-hybridization filter for viral bait set
        redundbaits      Run filter to remove highly similar (sequence space) baits
        humanwgc         Run filter to remove non-human virus targetting baits (uses Virus-Host DB)
        normalize        Generate a randomized subset bait set with a maximum number of sequences per virus family (EXPERIMENTAL)

        ---------------------------------------------------------------------------------------
        Statistics
        ---------------------------------------------------------------------------------------
        stats            Run all modules for bait statistics given a bait set FASTA file using default settings (WIP)
        proportion       Breakdown bait set by virus family and Balitmore classifers (i.e., +ssRNA, dsDNA, etc.)
        physprop         Calculate GC Content and melting temperature for a given FASTA file with baits
        coverage         Calculate depth and breadth of coverage for a given bait set FASTA (optionally provided reference)
        plots            Generate summary virus family coverage plots

        ---------------------------------------------------------------------------------------
        General
        ---------------------------------------------------------------------------------------
        help             Display this help page
        validate         Ensure all required packages and software are installed (will run through all functions with test data) (WIP)

```

