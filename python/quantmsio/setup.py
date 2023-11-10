import codecs
import os.path

from setuptools import find_packages
from setuptools import setup

with open("README.md", encoding="UTF-8") as fh:
    long_description = fh.read()


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), "r") as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    raise RuntimeError("Unable to find version string.")


setup(
    name="quantmsio",
    version=get_version("quantms_io/__init__.py"),
    author="BigBio Team",
    author_email="ypriverol@gmail.com",
    description="Convert multiple files into quantms.io format",
    long_description_content_type="text/markdown",
    long_description=long_description,
    license="'Apache 2.0",
    data_files=[("", ["LICENSE"])],
    package_data={"": ["*.xml"],},
    url="https://github.com/bigbio/quantms.io",
    packages=find_packages(),
    install_requires=["click", "pandas", "pandas_schema", "requests", "pytest", "pyyaml", "pyarrow",
                      "diskcache", "pyopenms", "scipy", "biopython", "seaborn"],
    entry_points={"console_scripts": ["quantmsio_cli = quantms_io.quantmsio_cli:quantms_io_main"]},
    platforms=["any"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="python multiomics proteomics quantms reanalysis",
    include_package_data=True,
    python_requires=">=3.6",
)