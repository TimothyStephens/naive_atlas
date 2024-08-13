from setuptools import setup
import versioneer  # script in directory

__author__ = "atlas: Silas Kieser & Joe Brown; naive_atlas: Timothy Stephens"
__copyright__ = "Copyright 2021, Silas Kieser"
__email__ = "atlas: silas.kieser@gmail.com & brwnjm@gmail.com; naive_atlas: ts942@sebs.rutgers.edu"
__license__ = "BSD-3"

# read the contents of your README file
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


setup(
    name="metagenome-naive_atlas",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    url="https://github.com/TimothyStephens/naive_atlas",
    license=__license__,
    author=__author__,
    author_email=__email__,
    zip_safe=False,
    description="NAIVE ATLAS - workflows for assembly, annotation, and genomic binning of mixed domain (prok, euk, virus) metagenomic and metatranscriptomic data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=["naive_atlas", "naive_atlas.init"],
    package_data={
        "": [
            "workflow",
        ]
    },
    data_files=[(".", ["README.md", "LICENSE.txt"])],
    include_package_data=True,
    install_requires=[],
    # install via conda: click, pandas, pyyaml, snakemake
    entry_points={"console_scripts": ["naive_atlas = naive_atlas.naive_atlas:cli"]},
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
)
