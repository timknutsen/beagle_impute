# rules/common.smk — helpers shared by Snakefile, Snakefile_accuracy and
# Snakefile_refpanel. Included first by each, so everything below is defined
# before any rule file needs it.

import functools
import os
import subprocess

import pandas as pd


@functools.lru_cache(maxsize=None)
def sorted_chromosomes(bfile):
    """Chromosome codes in a .bim, numeric codes first in numeric order.

    PLINK writes .bim tab-separated, but files that have passed through awk
    are often space-separated; \\s+ accepts both. X, Y, MT and scaffolds sort
    after the numeric codes rather than raising on int().
    """
    bim = pd.read_csv(f"{bfile}.bim", sep=r"\s+", header=None, usecols=[0], dtype=str)
    chroms = [c for c in bim[0].unique() if c not in ("0", "")]
    return sorted(chroms, key=lambda c: (0, int(c), "") if c.isdigit() else (1, 0, c))


def nested_config(section, key, default=None):
    """Resolve section.key from its flat --config alias first, then the YAML block.

    Snakemake rejects dotted keys on --config, so a nested setting is overridden
    from the CLI as section_key=value. The alias must win: the YAML files define
    every key, so reading the block first accepted a CLI override and silently
    ignored it.
    """
    if f"{section}_{key}" in config:
        return config[f"{section}_{key}"]
    if isinstance(config.get(section), dict) and key in config[section]:
        return config[section][key]
    return default


# -Xmx bounds the heap only; metaspace, thread stacks, GC structures and native
# buffers live outside it. Handing Java the full cgroup limit invites an OOM
# kill once the heap fills, so leave headroom.
_JAVA_HEAP_FRACTION = 0.85


def java_heap_mb(mem_mb):
    return max(1024, int(mem_mb * _JAVA_HEAP_FRACTION))


def beagle_map_param():
    """map=<file> for Beagle, or "" when beagle_params.map is unset.

    Without a map Beagle assumes 1 cM per Mb, so window/overlap are physical
    lengths and phasing ignores where the species recombines. Salmon
    recombination is concentrated near chromosome ends, which is where
    phase-switch errors cluster without one. The file is a PLINK .map
    (chrom, id, cM, bp; see scripts/linkage_map_to_plink.py); its chromosome
    codes must match the VCF CHROM column.
    """
    path = nested_config("beagle_params", "map", "")
    if not path:
        return ""
    if not os.path.isfile(path):
        raise ValueError(f"beagle_params.map={path!r} does not exist")
    return f"map={path}"


# JARs configured under bin/ are treated as auto-managed and fetched on first
# run. Point the config at any other path to use an existing JAR instead.
BEAGLE_URL  = "https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar"
BREF3_URL   = "https://faculty.washington.edu/browning/beagle/bref3.27Feb25.75f.jar"
CONFORM_URL = "https://faculty.washington.edu/browning/conform-gt/conform-gt.24May16.cee.jar"


def auto_download(path, url):
    if str(path).startswith("bin/") and not os.path.exists(path):
        os.makedirs("bin", exist_ok=True)
        print(f"Downloading {url} -> {path}")
        subprocess.run(["wget", "-q", "-O", path, url], check=True)


# Always run plink2 with --dog so non-human chromosome codes (1-38) are
# accepted: salmon (29), trout (32) and livestock fit, and human data still
# works. Prepended once here so every rule inherits it via plink_extra_flags.
_plink_flags = str(config.get("plink_extra_flags") or "")
if "--dog" not in _plink_flags.split():
    config["plink_extra_flags"] = ("--dog " + _plink_flags).strip()
