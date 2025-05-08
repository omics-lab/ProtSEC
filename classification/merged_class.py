

sub_to_super = {
    "Cell.membrane-M": "Cell.membrane",
    "Cytoplasm-Nucleus-U": "Cytoplasm",
    "Cytoplasm-S": "Cytoplasm",
    "Endoplasmic.reticulum-M": "Endoplasmic.reticulum",
    "Endoplasmic.reticulum-S": "Endoplasmic.reticulum",
    "Endoplasmic.reticulum-U": "Endoplasmic.reticulum",
    "Extracellular-S": "Extracellular",
    "Golgi.apparatus-M": "Golgi.apparatus",
    "Golgi.apparatus-S": "Golgi.apparatus",
    "Golgi.apparatus-U": "Golgi.apparatus",
    "Lysosome/Vacuole-M": "Lysosome/Vacuole",
    "Lysosome/Vacuole-S": "Lysosome/Vacuole",
    "Lysosome/Vacuole-U": "Lysosome/Vacuole",
    "Mitochondrion-M": "Mitochondrion",
    "Mitochondrion-S": "Mitochondrion",
    "Mitochondrion-U": "Mitochondrion",
    "Nucleus-M": "Nucleus",
    "Nucleus-S": "Nucleus",
    "Nucleus-U": "Nucleus",
    "Peroxisome-M": "Peroxisome",
    "Peroxisome-S": "Peroxisome",
    "Peroxisome-U": "Peroxisome",
    "Plastid-M": "Plastid",
    "Plastid-S": "Plastid",
    "Plastid-U": "Plastid", 
}


if __name__ == "__main__":
    s = {value for key, value in sub_t_super.items()}
    print(s, len(s))
