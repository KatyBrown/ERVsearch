import setuptools
import re

with open("README.md", "r") as fh:
    long_description = fh.read()

VERSIONFILE="ERVsearch/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

setuptools.setup(
     name='ervsearch',
     version=verstr,
     author="Katy Brown",
     author_email="kab84@cam.ac.uk",
     description="Tool for identifying endogenous retrovirus like regions in a set of sequences",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/KatyBrown/ERVsearch",
     packages=setuptools.find_packages(),
     install_requires=['numpy', 'matplotlib', 'ruffus', 'pandas', 'ete3'],
     include_package_data=True,
     scripts=['ERVsearch/ERVsearch'],
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent"]
 )
