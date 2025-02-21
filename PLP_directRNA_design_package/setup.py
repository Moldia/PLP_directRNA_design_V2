from setuptools import setup, find_packages

setup(
    name='plp_directrna_design',
    version='1.0.0',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'biopython'
    ],
    entry_points={
        "console_scripts": [
            "extract_features=PLP_directRNA_design.extract_features:main",
        ],
    },
    author='Nima Rafati, Marco Grillo, Victoria Muiliadi, Maria Escriva Conde, and Sergio Marco Salas',
    description='A package for designing probes for direct RNA sequencing.'
)