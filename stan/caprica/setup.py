import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

with open("caprica/_version.py", "r") as f:
    line = f.read().strip().replace(" ", "")
    version = line.split("=")[-1]

with open("requirements.txt", "r") as f:
    requirements = f.readlines()

setuptools.setup(
    name="caprica",
    version=version,
    author="Asher Wasserman",
    author_email="awasserman@xcures.com",
    description="GCTA simulation framework",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/adwasser/caprica",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
    ],
    install_requires=requirements,
    python_requires=">=3.6",
)
