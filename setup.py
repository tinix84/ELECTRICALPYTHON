# Import Necessary Files
import setuptools
import os

with open("README.md", "r") as fh:
	long_description = fh.read()
	
# Version Breakdown:
# MAJOR CHANGE . MINOR CHANGE

# Generate Setup Tools Argument
setuptools.setup(
	name="electricalpython",
	version="1.0",
	author="Joe Stanley",
	author_email="stan3926@vandals.uidaho.edu",
	description="Electrical Engineering Functions in Python",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/engineerjoe440/electrical-engineering-python",
	packages=setuptools.find_packages(),
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
	],
	project_urls={
        "Bug Tracker": "https://github.com/engineerjoe440/ELECTRICALPYTHON/issues",
        "Documentation": "https://engineerjoe440.github.io/ELECTRICALPYTHON/",
		}
)


# End of SETUP.PY