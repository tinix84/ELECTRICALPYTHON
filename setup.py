import setuptools

with open("README.md", "r") as fh:
	long_description = fh.read()

setuptools.setup(
	name="eepython",
	version="0.0.1",
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
)