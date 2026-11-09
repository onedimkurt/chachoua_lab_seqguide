"""SeqGuide: a guided sequencing analysis pipeline for non-expert users.

Core modules are pure Python so they can be tested without bioinformatics
tools installed. External tools (fastp, Unicycler, Flye, minimap2, samtools,
bcftools) are only called from `workflows.py`.
"""

__version__ = "0.4.0"
